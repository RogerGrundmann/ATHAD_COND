#pragma once

#include "cAtmosphereModel.h"
#include "Utils.h"

#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace AtomUtils;

class VelocityInitializer {
public:
    explicit VelocityInitializer(cAtmosphereModel& model)
        : m(model)
    {}

    void compute()
    {
        using namespace std;
        cout << endl << "      AGCM: init_velocities" << endl;
        printCellGeometry();

        const int n = nCells();

        // Radial u lives on the cell EDGES — the ascent/descent branches. Sign
        // alternates outward from the equator, so cell k is closed by a rising branch
        // on one side and a sinking branch on the other whatever n is.
        for(int k = 0; k <= n; k++){
            const double phi   = edgeLat(k);
            const double coeff = edgeRadialCoeff(k, phi);
            init_u(m.u, jN(phi), coeff);
            if(jS(phi) != jN(phi)) init_u(m.u, jS(phi), coeff);
        }

        // v and w live on the edges AND the cell cores. Both hemispheres take the SAME
        // coefficients — the southern sign flip happens in the fused pass below, which
        // is why nothing here is mirrored by hand. The southern index is derived from
        // the northern one rather than rounded independently, so invariant 1 survives
        // an n whose anchors do not land on integer indices (n = 4 puts an edge at
        // 22.5 degrees).
        for(int k = 0; k <= n; k++) setAnchor(edgeLat(k),   edgeAmp(k, n));
        for(int k = 0; k <  n; k++) setAnchor(centreLat(k), centreAmp(k, n));

        // Linear fill between consecutive anchors of each field. The inherited code
        // wrote these 24 calls out by hand in a fixed order; consecutive-pair chaining
        // reproduces exactly that set at n = 3, and is the only form that survives a
        // change of n.
        formChain(m.u, anchorList(false, false));
        formChain(m.v, anchorList(true,  false));
        formChain(m.w, anchorList(true,  true));

        // Zero land cells; non-dimensionalise air cells — single fused pass
        const double inv_u_0 = 1.0 / m.u_0;

        #pragma omp parallel for collapse(2) schedule(static)
        for (int i = 0; i < m.im; i++) {
            for (int k = 0; k < m.km; k++) {
                for (int j = 0; j < m.jm; j++) {
                    if (is_land(m.h, i, j, k)) {
                        m.u.x[i][j][k] = 0.0;
                        m.v.x[i][j][k] = 0.0;
                        m.w.x[i][j][k] = 0.0;
                    } else {
                        m.u.x[i][j][k] *= inv_u_0;
                        m.w.x[i][j][k] *= inv_u_0;
                        if (!m.use_NASA_velocity && j > 90) {
                            m.v.x[i][j][k] = -m.v.x[i][j][k] * inv_u_0;
                        } else {
                            m.v.x[i][j][k] *= inv_u_0;
                        }
                    }
                }
            }
        }
/*
        // Surface taper on v and w (the two HORIZONTAL components in this model's
        // (r,θ,φ) convention: v = meridional, w = zonal): linearly damp from the
        // local value at i=5 down to zero at i=0, so the lowest five layers carry no
        // horizontal wind at the ground reference and grow smoothly into the
        // prescribed profile above. u is left alone because u is the RADIAL/VERTICAL
        // velocity here (NOT the zonal jet — that is w); init_u already gives it a
        // small profile that ramps to zero at the surface.
        #pragma omp parallel for collapse(2) schedule(static)
        for (int j = 0; j < m.jm; j++) {
            for (int k = 0; k < m.km; k++) {
                for (int i = 0; i <= 5; i++) {
                    const double factor = static_cast<double>(i) / 5.0;
                    m.v.x[i][j][k] *= factor;
                    m.w.x[i][j][k] *= factor;
                }
            }
        }
*/
    cout << "      AGCM: init_velocities ended" << endl;
    }

    // Linear blend of u/v/w across j in [lat-3, lat+3].
    // Currently not called by compute() (dead code in the original),
    // but kept here as it logically belongs with velocity initialisation.
    void smooth_transition(int lat)
    {
        const int    start     = lat - 3;
        const int    end       = lat + 3;
        const double inv_range = 1.0 / (double)(end - start);

        #pragma omp parallel for collapse(2) schedule(static)
        for (int k = 0; k < m.km; k++) {
            for (int i = 0; i < m.im; i++) {
                const double u_start = m.u.x[i][start][k];
                const double v_start = m.v.x[i][start][k];
                const double w_start = m.w.x[i][start][k];
                const double u_slope = (m.u.x[i][end][k] - u_start) * inv_range;
                const double v_slope = (m.v.x[i][end][k] - v_start) * inv_range;
                const double w_slope = (m.w.x[i][end][k] - w_start) * inv_range;
                for (int j = start; j <= end; j++) {
                    const double t     = (double)(j - start);
                    m.u.x[i][j][k] = u_slope * t + u_start;
                    m.v.x[i][j][k] = v_slope * t + v_start;
                    m.w.x[i][j][k] = w_slope * t + w_start;
                }
            }
        }
    }

private:
    cAtmosphereModel& m;


    // cell_lat_scale — compress the prescribed cell latitudes toward the equator.
    //
    // Every anchor below is written at an EARTH latitude: Hadley 15 deg, Ferrel 45, polar
    // 75, with the trade/westerly nodes between them. Those latitudes are a consequence of
    // Earth's thermal Rossby number, and this atmosphere's is 12.5x smaller —
    // Ro_T = g*H*(dtheta/theta)/(omega^2 a^2) = 0.0048 against 0.0598, because rotation is
    // 4.35x faster (omega^2 18.9x) and the FRACTIONAL equator-pole contrast is 4.7x weaker
    // (50 K on 1500 K against 45 K on 288 K), only partly offset by a 7x deeper atmosphere.
    // Held-Hou then puts the direct cell's edge near 5 deg here against 18 deg for Earth,
    // i.e. cells roughly a third as wide — hence the 0.33 default. A hot surface is not a
    // strongly DIFFERENTIALLY heated one, which is why the higher energy content narrows
    // the circulation instead of widening it.
    //
    // Set 1.0 to restore Earth's latitudes, which is what every run before README item 31
    // used. s now scales the HADLEY EDGE ONLY — see edgeLat below for why scaling every
    // anchor uniformly does not tile the hemisphere. The measured cell decay and the
    // derivation are in param.py.
    double latScale() const { return m.cell_lat_scale; }

    // ---- cell geometry ------------------------------------------------------
    // n cells per hemisphere. Cell k (k = 0 at the equator) runs from edge k to
    // edge k+1 with its core at the midpoint of the two.
    //
    // THE CELLS MUST TILE THE HEMISPHERE. Scaling every anchor by s and letting js()
    // pin the pole does not: at s = 0.33 it gives edges 0/10/20/90, so the three cells
    // are 10, 10 and *70* degrees wide and the whole extratropics is one linear ramp
    // between the 20 deg anchor and the pole. Adding cells made it worse, not better --
    // n = 4 packed four narrow cells into 0-22 deg and left a 68 deg cell behind. That
    // was the committed behaviour from item 31 onward, and it is why item 31's scan
    // measured a tropical cell against an extratropics that had been deleted rather
    // than narrowed.
    //
    // Held-Hou constrains the DIRECT cell's edge and says nothing about the
    // extratropical bands, which the Rhines scale sets independently. So s scales the
    // Hadley edge only, and the remaining band is tiled by the other n-1 cells:
    //
    //     edge(0) = 0
    //     edge(1) = 30 * s                                    <- Earth's Hadley edge, scaled
    //     edge(k) = edge(1) + (90 - edge(1))*(k-1)/(n-1)
    //     core(k) = midpoint of edge(k), edge(k+1)
    //
    // At s = 1, n = 3 this is 0/30/60/90 with cores 15/45/75 -- Earth exactly, midpoints
    // and all. At s = 0.33 the extratropical bands come out 26.7 deg wide at n = 4 and
    // 20 deg at n = 5, against the Rhines scale's ~22 deg, which is the first time the
    // Held-Hou and Rhines estimates have agreed on a layout.
    //
    // jN must NOT go through js() any more: the scaling is already in edgeLat, and js
    // would apply it a second time. The pole needs no pinning either, because edge(n) is
    // exactly 90 by construction.
    int    nCells()          const { return m.n_cells_hemisphere; }
    double hadleyEdge()      const { return 30.0 * latScale(); }
    double edgeLat(int k)    const {
        if(k <= 0)         return 0.0;
        if(nCells() <= 1)  return 90.0;
        const double h = hadleyEdge();
        return h + (90.0 - h) * (double)(k - 1) / (double)(nCells() - 1);
    }
    double centreLat(int k)  const { return 0.5 * (edgeLat(k) + edgeLat(k + 1)); }
    int    jN(double phi)    const {
        int j = (int)std::lround(90.0 - phi);
        if(j < 0)          j = 0;
        if(j > m.jm - 1)   j = m.jm - 1;
        return j;
    }
    int    jS(double phi)    const { return (m.jm - 1) - jN(phi); }

    // ---- prescribed amplitudes ----------------------------------------------
    // Earth's inherited table, re-expressed by ROLE rather than by latitude literal.
    // has_w is false for the polar core, which the inherited code deliberately left
    // to interpolation rather than anchoring.
    struct Amp { double v_trop, v_surf, w_trop, w_surf; bool has_w; };

    static Amp edgeEquator()  { return {  0.0,  0.0, -3.0, -5.0, true  }; }
    static Amp edgeSubtrop()  { return {  0.0,  0.5, 30.0, -1.0, true  }; }
    static Amp edgeSubpolar() { return { -0.2,  0.0, 10.0,  6.0, true  }; }
    static Amp edgePole()     { return {  0.5,  0.0,  0.0,  0.0, true  }; }
    static Amp cellHadley()   { return { -3.0,  3.5,  5.0, -7.0, true  }; }
    static Amp cellFerrel()   { return {  4.0, -1.5, 15.0, 10.0, true  }; }
    static Amp cellPolar()    { return {  0.5,  0.6,  0.0,  0.0, false }; }

    // Cell 0 keeps the Hadley template, cell n-1 the polar one, and everything
    // between is a Ferrel copy — so raising n INSERTS bands in mid-latitudes, where
    // the Rhines argument says the deformation scale shrinks, instead of at the pole.
    // At n = 3 this is Hadley / Ferrel / polar, unchanged. See param.py for why the
    // direct-indirect alternation cannot be used as the rule instead: the inherited
    // polar cell carries the Ferrel's sense, not the Hadley's.
    static Amp edgeAmp(int k, int n){
        if(k == 0) return edgeEquator();
        if(k == n) return edgePole();
        if(k == 1) return edgeSubtrop();
        return edgeSubpolar();
    }
    static Amp centreAmp(int k, int n){
        if(k == 0)     return cellHadley();
        if(k == n - 1) return cellPolar();
        return cellFerrel();
    }

    // The four inherited radial magnitudes are 0.02894*(1 - phi/150) to every figure
    // they are written with, but they are kept as literals so n = 3 stays bit-identical;
    // the formula only supplies edges Earth has no value for.
    static double edgeRadialCoeff(int k, double phi){
        double mag;
        if      (phi ==  0.0) mag = 0.02894;
        else if (phi == 30.0) mag = 0.02315;
        else if (phi == 60.0) mag = 0.01736;
        else if (phi == 90.0) mag = 0.011574;
        else                  mag = 0.02894 * (1.0 - phi / 150.0);
        return (k % 2 == 0) ? mag : -mag;
    }

    void setAnchor(double phi, const Amp& a){
        const int jn = jN(phi), jsouth = jS(phi);
        init_v_or_w(m.v, jn, a.v_trop, a.v_surf);
        if(jsouth != jn) init_v_or_w(m.v, jsouth, a.v_trop, a.v_surf);
        if(a.has_w){
            init_v_or_w(m.w, jn, a.w_trop, a.w_surf);
            if(jsouth != jn) init_v_or_w(m.w, jsouth, a.w_trop, a.w_surf);
        }
    }

    // Sorted, de-duplicated anchor indices for one field, pole to pole.
    std::vector<int> anchorList(bool with_centres, bool skip_polar_centre) const {
        const int n = nCells();
        std::vector<double> phis;
        for(int k = 0; k <= n; k++) phis.push_back(edgeLat(k));
        if(with_centres)
            for(int k = 0; k < n; k++){
                if(skip_polar_centre && k == n - 1) continue;
                phis.push_back(centreLat(k));
            }
        std::vector<int> j;
        for(double p : phis){
            j.push_back(jN(p));
            if(jS(p) != jN(p)) j.push_back(jS(p));
        }
        std::sort(j.begin(), j.end());
        j.erase(std::unique(j.begin(), j.end()), j.end());
        return j;
    }

    void formChain(Array& a, const std::vector<int>& anchors){
        for(size_t i = 0; i + 1 < anchors.size(); i++)
            if(anchors[i + 1] > anchors[i]) form_diagonals(a, anchors[i], anchors[i + 1]);
    }

    void printCellGeometry() const {
        using namespace std;
        const int n = nCells();
        cout << "      n_cells_hemisphere = " << n << ", cell_lat_scale = " << latScale()
             << (n == 3 && latScale() == 1.0 ? "   <- Earth's 3-cell layout" : "") << endl;
        cout << "      cell edges at ";
        for(int k = 0; k <= n; k++)
            cout << (90 - jN(edgeLat(k))) << (k < n ? " / " : "");
        cout << " deg N; cores at ";
        for(int k = 0; k < n; k++)
            cout << (90 - jN(centreLat(k))) << (k + 1 < n ? " / " : "");
        cout << " deg; widths ";
        for(int k = 0; k < n; k++)
            cout << (jN(edgeLat(k)) - jN(edgeLat(k + 1))) << (k + 1 < n ? " / " : "");
        cout << " deg" << endl;
        cout << "      cell_amp_mode  = " << m.cell_amp_mode << ": v x " << ampScale(m.v)
             << ", w x " << ampScale(m.w) << ", radial x 1 (continuity)"
             << (m.cell_amp_mode == 0 && latScale() != 1.0
                 ? "   <- amplitudes unscaled, meridional shear is 1/s x Earth's" : "")
             << endl;
    }

    // cell_amp_mode — scale the prescribed AMPLITUDES with the latitude scale.
    //
    // latScale() moves the anchors and nothing else, so at s = 0.33 the same velocity
    // change is interpolated across a third of the latitude span and every meridional
    // gradient in the initial state is 3x Earth's. These modes put the amplitudes on the
    // same footing as the anchors; all of them are a no-op at s = 1. See param.py for the
    // continuity and angular-momentum derivations.
    //
    // The RADIAL amplitudes (ua_00 .. ua_90 in init_u) are deliberately never scaled:
    // continuity makes du_r/dz invariant when v and the latitude scale shrink together.
    double ampScale(const Array& a) const {
        const double s = latScale();
        switch(m.cell_amp_mode){
            case 1: return s;                                 // continuity, both components
            case 2: return (&a == &m.v) ? s : s * s;          // continuity / angular momentum
            default: return 1.0;                              // off
        }
    }
    int js(int j) const {
        if(j <= 0) return 0;
        if(j >= m.jm - 1) return m.jm - 1;
        const int jeq = (m.jm - 1) / 2;
        return (int)std::lround(jeq + (double)(j - jeq) * latScale());
    }

    // j_earth is the UNSCALED anchor index — the case labels below are Earth's grid
    // indices, so the lookup has to happen before js() moves the anchor. Passing the
    // scaled index in (as compute() used to) makes the switch miss: at
    // cell_lat_scale = 0.33 the anchors land on 60/70/80/90/100/110/120, of which only
    // 90 is a real case, 60 and 120 pick up a NEIGHBOUR's coefficient (the poles got
    // +ua_60 instead of -ua_90 — sign flipped and 1.5x too large) and the other four
    // fall through `default: return` and are silently never written. At 0.75 six of the
    // seven miss. form_diagonals then interpolates across anchors that were never set.
    //
    // The write location and the tropopause layer must still use the SCALED index,
    // which is the whole point of the knob: same coefficient, moved in latitude.
    // j is the SCALED write index; coeff comes from edgeRadialCoeff, which is keyed to
    // the cell edge rather than to a grid index. The switch(j) this replaces took Earth's
    // grid indices as case labels while compute() passed the scaled index, so at
    // cell_lat_scale != 1 most anchors missed their case entirely — see eaf6e8d.
    void init_u(Array& u, int j, double coeff)
    {
        const int    tl           = m.get_tropopause_layer(j);
        const double tropo_h      = m.get_layer_height(tl);
        const double half_tropo_h = tropo_h / 3.0;
        const double inv_ascent   = 3.0 / half_tropo_h;
        const double inv_descent  = 1.0 / half_tropo_h;

        #pragma omp parallel for schedule(static)
        for (int k = 0; k < m.km; k++) {
            for (int i = 0; i < tl; i++) {
                const double h     = m.get_layer_height(i);
                const double ratio = (h < half_tropo_h)
                    ? h * inv_ascent
                    : (tropo_h - h) * inv_descent;
                u.x[i][j][k] = coeff * ratio;
            }
        }
    }

    void init_v_or_w(Array& v_or_w, int j, double coeff_trop, double coeff_sl)
    {
        const int    tl          = m.get_tropopause_layer(j);
        const double inv_tropo_h = 1.0 / m.get_layer_height(tl);

        // Amplitudes scale with the anchors (cell_amp_mode); the NASA-surface branch below
        // reads an observed value out of the field instead, so it must not be rescaled.
        const double amp = ampScale(v_or_w);
        coeff_trop *= amp;
        coeff_sl   *= amp;

        #pragma omp parallel for schedule(static)
        for (int k = 0; k < m.km; k++) {
            double sl = coeff_sl;
            if (m.use_NASA_velocity && is_ocean_surface(m.h, 0, j, k)) {
                sl = v_or_w.x[0][j][k];
            }
            const double slope = (coeff_trop - sl) * inv_tropo_h;

            for (int i = 0; i < tl; i++) {
                v_or_w.x[i][j][k] = slope * m.get_layer_height(i) + sl;
            }
        }

        init_v_or_w_above_tropopause(v_or_w, j, coeff_trop);
    }

    void init_v_or_w_above_tropopause(Array& v_or_w, int j, double coeff)
    {
        const int tl = m.get_tropopause_layer(j);
        if (tl >= m.im - 1) return;

        const double h_top     = m.get_layer_height(m.im - 1);
        const double inv_range = 1.0 / (h_top - m.get_layer_height(tl));

        #pragma omp parallel for schedule(static)
        for (int k = 0; k < m.km; k++) {
            for (int i = tl; i < m.im; i++) {
                v_or_w.x[i][j][k] = coeff * (h_top - m.get_layer_height(i)) * inv_range;
            }
        }
    }

    void form_diagonals(Array& a, int start, int end)
    {
        const double inv_range = 1.0 / (double)(end - start);

        #pragma omp parallel for collapse(2) schedule(static)
        for (int k = 0; k < m.km; k++) {
            for (int i = 0; i < m.im; i++) {
                const double a_start = a.x[i][start][k];
                const double slope   = (a.x[i][end][k] - a_start) * inv_range;
                for (int j = start; j < end; j++) {
                    a.x[i][j][k] = slope * (double)(j - start) + a_start;
                }
            }
        }
    }
};
