#pragma once

#include "cAtmosphereModel.h"
#include "Utils.h"

#include <iostream>
#include <algorithm>
#include <cmath>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace AtomUtils;

class BC_Atm {
public:
    explicit BC_Atm(cAtmosphereModel& model)
        : m(model)
    {}

    // ------------------------------------------------------------------
    void bcSolidGround()
    {
        using namespace std;
        cout << endl << endl << endl << "      AGCM: BC_SolidGround" << endl;

        // Reference upward long-wave flux for reset land cells — the radiative
        // counterpart of the t = 1.0 reference (T = t_0), i.e. σ·t_0⁴ [W/m²].
        // radiation.x = σ·(t·t_0)⁴ everywhere (see MultiLayerRadiation), so a solid
        // cell whose t is reset to 1.0 must have its radiation reset to σ·t_0⁴ to
        // stay consistent and not leave a stale flux inside the topography.
        const double rad_ref = m.sigma * std::pow(m.t_0, 4.0);

        // Multi-direction averaging from face-neighbouring air cells.
        //
        // Replaces the previous (4/3)·x_near − (1/3)·x_far one-sided Neumann
        // extrapolation, which amplified any boundary noise by 4/3 and overshot
        // at steep topography.  Observed problem: the Himalaya peak, where a single
        // land ghost cell has air neighbours in multiple directions of very different
        // air-column heights, drove spurious p_dyn and velocity anomalies that
        // propagated downstream and destabilised the whole NH at iter ~120.
        //
        // The new method: for each land ghost cell, look at the 6 face neighbours
        // (i±1, j±1, k±1 with k periodic).  Average the values from whichever of
        // those neighbours are air cells.  No amplification (factor 1.0 by
        // construction), direction-independent (a peak with air on top, on the east
        // face, and on the south face gets all three contributions equally weighted),
        // and naturally degrades to plain copy when only one neighbour is air.
        //
        // Non-negative fields (tke, dis, nue, c, cloud, ice, gr, q_*, E_*, D_*, M_*)
        // are clamped to [0, ∞); the averaging of non-negative values can only return
        // non-negative results, so this clamp is defensive.
        //
        // p_dyn was previously held at 0 here as a Dirichlet pin.  It now averages too,
        // so the pressure varies smoothly across the topographic surface.  The Poisson
        // anchor is provided by (a) the still-zero p_dyn at fully-buried interior cells
        // (Pass 1) and (b) the PressureSolverAtm i=0-over-land pin.
        auto average_from_air_neighbors = [&](int i, int j, int k) {
            int kp = (k + 1)        % m.km;
            int km = (k - 1 + m.km) % m.km;
            bool a_ip = (i + 1 < m.im) && is_air(m.h, i+1, j, k);
            bool a_im = (i - 1 >= 0)    && is_air(m.h, i-1, j, k);
            bool a_jp = (j + 1 < m.jm) && is_air(m.h, i, j+1, k);
            bool a_jm = (j - 1 >= 0)    && is_air(m.h, i, j-1, k);
            bool a_kp = is_air(m.h, i, j, kp);
            bool a_km = is_air(m.h, i, j, km);

            int count = (int)a_ip + (int)a_im + (int)a_jp + (int)a_jm + (int)a_kp + (int)a_km;
            if (count == 0) return;
            const double inv = 1.0 / static_cast<double>(count);

            auto avg = [&](const Array& f) -> double {
                double sum = 0.0;
                if (a_ip) sum += f.x[i+1][j][k];
                if (a_im) sum += f.x[i-1][j][k];
                if (a_jp) sum += f.x[i][j+1][k];
                if (a_jm) sum += f.x[i][j-1][k];
                if (a_kp) sum += f.x[i][j][kp];
                if (a_km) sum += f.x[i][j][km];
                return sum * inv;
            };

            m.tke.x[i][j][k] = std::max(0.0, avg(m.tke));
            m.dis.x[i][j][k] = std::max(0.0, avg(m.dis));
            m.nue.x[i][j][k] = std::max(0.0, avg(m.nue));

            m.t.x[i][j][k]   = avg(m.t);
            m.radiation.x[i][j][k] = avg(m.radiation);

            m.p_dyn.x[i][j][k] = avg(m.p_dyn);

            m.r_dry.x[i][j][k]   = avg(m.r_dry);
            m.r_humid.x[i][j][k] = avg(m.r_humid);

            m.c.x[i][j][k]     = std::max(0.0, avg(m.c));
            m.cloud.x[i][j][k] = std::max(0.0, avg(m.cloud));
            m.ice.x[i][j][k]   = std::max(0.0, avg(m.ice));
            m.gr.x[i][j][k]    = std::max(0.0, avg(m.gr));

            m.co2.x[i][j][k]   = avg(m.co2);

            m.S_c_c.x[i][j][k] = avg(m.S_c_c);
            m.S_v.x[i][j][k]   = avg(m.S_v);
            m.S_c.x[i][j][k]   = avg(m.S_c);
            m.S_i.x[i][j][k]   = avg(m.S_i);
            m.S_r.x[i][j][k]   = avg(m.S_r);
            m.S_s.x[i][j][k]   = avg(m.S_s);

            m.q_v_u.x[i][j][k] = std::max(0.0, avg(m.q_v_u));
            m.q_c_u.x[i][j][k] = std::max(0.0, avg(m.q_c_u));
            m.q_v_d.x[i][j][k] = std::max(0.0, avg(m.q_v_d));

            m.E_u.x[i][j][k] = std::max(0.0, avg(m.E_u));
            m.E_d.x[i][j][k] = std::max(0.0, avg(m.E_d));
            m.D_u.x[i][j][k] = std::max(0.0, avg(m.D_u));
            m.D_d.x[i][j][k] = std::max(0.0, avg(m.D_d));
            m.M_u.x[i][j][k] = std::max(0.0, avg(m.M_u));
            m.M_d.x[i][j][k] = std::max(0.0, avg(m.M_d));
        };

        // Pass 1 — apply solid-wall conditions to every land cell from i=0 up to
        // and including the topographic surface layer i_mount.
        // Two distinct cases are handled in separate branches:
        //   interior : i < i_mount AND no air neighbour in j or k
        //              (fully buried — cannot affect the air column in any direction)
        //   surface  : i == i_mount OR adjacent to air in j or k
        //              (topographic top or cliff face — boundary to the air column)
        #pragma omp parallel for schedule(static)
        for (int j = 0; j < m.jm; j++) {
            for (int k = 0; k < m.km; k++) {
                int i_mount = m.i_topography[j][k];

                for (int i = 0; i <= i_mount; i++) {
                    if (!is_land(m.h, i, j, k)) continue;

                    // Check for air neighbours in j (latitude) and k (longitude).
                    // k is periodic: the Greenwich seam (k=0) wraps to k=km-1.
                    bool adj_air_j =
                        (j > 0       && is_air(m.h, i, j-1, k)) ||
                        (j < m.jm-1  && is_air(m.h, i, j+1, k));
                    int k_prev = (k > 0) ? k-1 : m.km-1;
                    int k_next = (k < m.km-1) ? k+1 : 0;
                    bool adj_air_k = is_air(m.h, i, j, k_prev) || is_air(m.h, i, j, k_next);

                    bool on_surface = (i == i_mount) || adj_air_j || adj_air_k;

                    if (!on_surface) {
                        // Fully buried interior cell — zero every prognostic field so
                        // the orography cannot contaminate the overlying air column.

                        // No-slip velocity (current and next time level).
                        m.u.x[i][j][k]  = 0.0;
                        m.v.x[i][j][k]  = 0.0;
                        m.w.x[i][j][k]  = 0.0;
                        m.un.x[i][j][k] = 0.0;
                        m.vn.x[i][j][k] = 0.0;
                        m.wn.x[i][j][k] = 0.0;

                        m.t.x[i][j][k]     = 1.0;
                        // Reset the time-level-n temperature too. The RK4 loop carries no
                        // land mask: it integrates every cell i=1..im-1 from tn, so it
                        // overwrites the t.x=1.0 reset above and reads tn instead. Without
                        // resetting tn, solid cliff-face cells (adjacent to ocean in j/k)
                        // accumulate an unbounded cold anomaly in tn each step (observed as a
                        // growing sub-terrain cold pool at the steep BC/Alaska coast). Mirrors
                        // the un/vn/wn and cn/cloudn/icen/grn resets elsewhere in this branch.
                        m.tn.x[i][j][k]    = 1.0;
                        m.radiation.x[i][j][k] = rad_ref;
                        m.p_dyn.x[i][j][k] = 0.0;

                        // Density reset to dry-air reference. CO2: reset buried cells to the
                        // column's surface value co2.x[0] (not the co2_0 constant) so (a) no
                        // artificial sub-terrain->air jump packs contour lines at mountains, and
                        // (b) it stays consistent with co2_scale (co2.x[0] is already scaled).
                        m.r_dry.x[i][j][k]   = m.r_air;
                        m.r_humid.x[i][j][k] = m.r_air;
                        m.co2.x[i][j][k]     = m.co2.x[0][j][k];

                        // Microphysics (current and time-level n+1).
                        m.c.x[i][j][k]      = 0.0;
                        m.cloud.x[i][j][k]  = 0.0;
                        m.ice.x[i][j][k]    = 0.0;
                        m.gr.x[i][j][k]     = 0.0;
                        m.cn.x[i][j][k]     = 0.0;
                        m.cloudn.x[i][j][k] = 0.0;
                        m.icen.x[i][j][k]   = 0.0;
                        m.grn.x[i][j][k]    = 0.0;

                        // Microphysical source/sink terms.
                        m.S_c_c.x[i][j][k] = 0.0;
                        m.S_v.x[i][j][k]   = 0.0;
                        m.S_c.x[i][j][k]   = 0.0;
                        m.S_i.x[i][j][k]   = 0.0;
                        m.S_r.x[i][j][k]   = 0.0;
                        m.S_s.x[i][j][k]   = 0.0;

                        // All body forces zero inside the solid body.
                        m.BuoyancyForce.x[i][j][k]    = 0.0;
                        m.PresGradForce.x[i][j][k]    = 0.0;
                        m.CoriolisForce.x[i][j][k]    = 0.0;
                        m.CentrifugalForce.x[i][j][k] = 0.0;

                        m.Q_Latent.x[i][j][k]   = 0.0;
                        m.Q_Sensible.x[i][j][k] = 0.0;

                        // Turbulence: zero all scalars inside solid cells.
                        m.tke.x[i][j][k]        = 0.0;
                        m.tken.x[i][j][k]       = 0.0;
                        m.dis.x[i][j][k]        = 0.0;
                        m.disn.x[i][j][k]       = 0.0;
                        m.nue.x[i][j][k]        = 0.0;
                        m.prod.x[i][j][k]       = 0.0;
                        m.tke_source.x[i][j][k] = 0.0;
                        m.dis_source.x[i][j][k] = 0.0;

                        // Moist-convection scalars: dry static energy set to reference (1),
                        // all convective velocities and moisture quantities zeroed.
                        m.s.x[i][j][k]   = 1.0;
                        m.s_u.x[i][j][k] = 1.0;
                        m.s_d.x[i][j][k] = 1.0;
                        m.u_u.x[i][j][k] = 0.0;
                        m.u_d.x[i][j][k] = 0.0;
                        m.v_u.x[i][j][k] = 0.0;
                        m.v_d.x[i][j][k] = 0.0;
                        m.w_u.x[i][j][k] = 0.0;
                        m.w_d.x[i][j][k] = 0.0;

                        // Convective condensate and precipitation conversion terms.
                        m.c_u.x[i][j][k] = 0.0;
                        m.e_d.x[i][j][k] = 0.0;
                        m.e_l.x[i][j][k] = 0.0;
                        m.e_p.x[i][j][k] = 0.0;
                        m.g_p.x[i][j][k] = 0.0;

                        m.q_c_u.x[i][j][k] = 0.0;
                        m.q_v_u.x[i][j][k] = 0.0;
                        m.q_v_d.x[i][j][k] = 0.0;

                        // Mass-flux and entrainment/detrainment rates.
                        m.D_d.x[i][j][k] = 0.0;
                        m.D_u.x[i][j][k] = 0.0;
                        m.E_d.x[i][j][k] = 0.0;
                        m.E_u.x[i][j][k] = 0.0;
                        m.M_d.x[i][j][k] = 0.0;
                        m.M_u.x[i][j][k] = 0.0;

                        // Moist-convection tendencies (temperature, humidity, momentum).
                        m.MC_q.x[i][j][k] = 0.0;
                        m.MC_t.x[i][j][k] = 0.0;
                        m.MC_v.x[i][j][k] = 0.0;
                        m.MC_w.x[i][j][k] = 0.0;
                    } else {
                        // Surface cell — either the topographic top (i == i_mount) or a
                        // cliff face (i < i_mount but adjacent to air in j or k).

                        if (!m.inviscid_phase) {
                            // No-slip velocity (current and next time level).
                            m.u.x[i][j][k]  = 0.0;
                            m.v.x[i][j][k]  = 0.0;
                            m.w.x[i][j][k]  = 0.0;
                            m.un.x[i][j][k] = 0.0;
                            m.vn.x[i][j][k] = 0.0;
                            m.wn.x[i][j][k] = 0.0;
                        } else {
                            // Free-slip (reflecting) wall: wall-normal component zero, tangential
                            // components mirrored from the adjacent air cell. Direction of "normal"
                            // is chosen by which neighbour is air. Priority: i+1 (top of column),
                            // then j-faces, then k-faces. Only one face is handled per cell — good
                            // enough for the spin-up phase; cliff corners get the i-face treatment.
                            bool air_above = (i + 1 <= m.im - 1) && is_air(m.h, i+1, j, k);
                            if (air_above) {
                                // Top of topographic column: normal is i. Zero u, mirror v,w.
                                m.u.x[i][j][k]  = 0.0;
                                m.v.x[i][j][k]  = m.v.x[i+1][j][k];
                                m.w.x[i][j][k]  = m.w.x[i+1][j][k];
                            } else if (j > 0 && is_air(m.h, i, j-1, k)) {
                                m.u.x[i][j][k]  = m.u.x[i][j-1][k];
                                m.v.x[i][j][k]  = 0.0;
                                m.w.x[i][j][k]  = m.w.x[i][j-1][k];
                            } else if (j < m.jm-1 && is_air(m.h, i, j+1, k)) {
                                m.u.x[i][j][k]  = m.u.x[i][j+1][k];
                                m.v.x[i][j][k]  = 0.0;
                                m.w.x[i][j][k]  = m.w.x[i][j+1][k];
                            } else {
                                int k_prev2 = (k > 0)        ? k - 1 : m.km - 1;
                                int k_next2 = (k < m.km - 1) ? k + 1 : 0;
                                if (is_air(m.h, i, j, k_prev2)) {
                                    m.u.x[i][j][k]  = m.u.x[i][j][k_prev2];
                                    m.v.x[i][j][k]  = m.v.x[i][j][k_prev2];
                                    m.w.x[i][j][k]  = 0.0;
                                } else if (is_air(m.h, i, j, k_next2)) {
                                    m.u.x[i][j][k]  = m.u.x[i][j][k_next2];
                                    m.v.x[i][j][k]  = m.v.x[i][j][k_next2];
                                    m.w.x[i][j][k]  = 0.0;
                                } else {
                                    // Cell flagged as surface but no air neighbour — should not
                                    // occur given on_surface check; fall back to no-slip.
                                    m.u.x[i][j][k]  = 0.0;
                                    m.v.x[i][j][k]  = 0.0;
                                    m.w.x[i][j][k]  = 0.0;
                                }
                            }
                            m.un.x[i][j][k] = m.u.x[i][j][k];
                            m.vn.x[i][j][k] = m.v.x[i][j][k];
                            m.wn.x[i][j][k] = m.w.x[i][j][k];
                        }

                        m.t.x[i][j][k]     = 1.0;
                        // Reset the time-level-n temperature too. The RK4 loop carries no
                        // land mask: it integrates every cell i=1..im-1 from tn, so it
                        // overwrites the t.x=1.0 reset above and reads tn instead. Without
                        // resetting tn, solid cliff-face cells (adjacent to ocean in j/k)
                        // accumulate an unbounded cold anomaly in tn each step (observed as a
                        // growing sub-terrain cold pool at the steep BC/Alaska coast). Mirrors
                        // the un/vn/wn and cn/cloudn/icen/grn resets elsewhere in this branch.
                        m.tn.x[i][j][k]    = 1.0;
                        m.radiation.x[i][j][k] = rad_ref;
                        m.p_dyn.x[i][j][k] = 0.0;

                        // Density reset to dry-air reference. CO2: reset buried cells to the
                        // column's surface value co2.x[0] (not the co2_0 constant) so (a) no
                        // artificial sub-terrain->air jump packs contour lines at mountains, and
                        // (b) it stays consistent with co2_scale (co2.x[0] is already scaled).
                        m.r_dry.x[i][j][k]   = m.r_air;
                        m.r_humid.x[i][j][k] = m.r_air;
                        m.co2.x[i][j][k]     = m.co2.x[0][j][k];

                        // Microphysics (current and time-level n+1).
                        m.c.x[i][j][k]      = 0.0;
                        m.cloud.x[i][j][k]  = 0.0;
                        m.ice.x[i][j][k]    = 0.0;
                        m.gr.x[i][j][k]     = 0.0;
                        m.cn.x[i][j][k]     = 0.0;
                        m.cloudn.x[i][j][k] = 0.0;
                        m.icen.x[i][j][k]   = 0.0;
                        m.grn.x[i][j][k]    = 0.0;

                        // Microphysical source/sink terms.
                        m.S_c_c.x[i][j][k] = 0.0;
                        m.S_v.x[i][j][k]   = 0.0;
                        m.S_c.x[i][j][k]   = 0.0;
                        m.S_i.x[i][j][k]   = 0.0;
                        m.S_r.x[i][j][k]   = 0.0;
                        m.S_s.x[i][j][k]   = 0.0;

                        // All body forces zero inside the solid body.
                        m.BuoyancyForce.x[i][j][k]    = 0.0;
                        m.PresGradForce.x[i][j][k]    = 0.0;
                        m.CoriolisForce.x[i][j][k]    = 0.0;
                        m.CentrifugalForce.x[i][j][k] = 0.0;

                        m.Q_Latent.x[i][j][k]   = 0.0;
                        m.Q_Sensible.x[i][j][k] = 0.0;

                        // Turbulence: zero all scalars at solid surfaces.
                        m.tke.x[i][j][k]        = 0.0;
                        m.tken.x[i][j][k]       = 0.0;
                        m.dis.x[i][j][k]        = 0.0;
                        m.disn.x[i][j][k]       = 0.0;
                        m.nue.x[i][j][k]        = 0.0;
                        m.prod.x[i][j][k]       = 0.0;
                        m.tke_source.x[i][j][k] = 0.0;
                        m.dis_source.x[i][j][k] = 0.0;

                        // Moist-convection scalars: dry static energy set to reference (1),
                        // all convective velocities and moisture quantities zeroed.
                        m.s.x[i][j][k]   = 1.0;
                        m.s_u.x[i][j][k] = 1.0;
                        m.s_d.x[i][j][k] = 1.0;
                        m.u_u.x[i][j][k] = 0.0;
                        m.u_d.x[i][j][k] = 0.0;
                        m.v_u.x[i][j][k] = 0.0;
                        m.v_d.x[i][j][k] = 0.0;
                        m.w_u.x[i][j][k] = 0.0;
                        m.w_d.x[i][j][k] = 0.0;

                        // Convective condensate and precipitation conversion terms.
                        m.c_u.x[i][j][k] = 0.0;
                        m.e_d.x[i][j][k] = 0.0;
                        m.e_l.x[i][j][k] = 0.0;
                        m.e_p.x[i][j][k] = 0.0;
                        m.g_p.x[i][j][k] = 0.0;

                        m.q_c_u.x[i][j][k] = 0.0;
                        m.q_v_u.x[i][j][k] = 0.0;
                        m.q_v_d.x[i][j][k] = 0.0;

                        // Mass-flux and entrainment/detrainment rates.
                        m.D_d.x[i][j][k] = 0.0;
                        m.D_u.x[i][j][k] = 0.0;
                        m.E_d.x[i][j][k] = 0.0;
                        m.E_u.x[i][j][k] = 0.0;
                        m.M_d.x[i][j][k] = 0.0;
                        m.M_u.x[i][j][k] = 0.0;

                        // Moist-convection tendencies (temperature, humidity, momentum).
                        m.MC_q.x[i][j][k] = 0.0;
                        m.MC_t.x[i][j][k] = 0.0;
                        m.MC_v.x[i][j][k] = 0.0;
                        m.MC_w.x[i][j][k] = 0.0;
                    }
                }
            }
        }

        // Pass 2 — at every land cell that touches the air column in any of its 6
        // face directions (i±1, j±1, k±1 with k periodic), replace its scalar fields
        // with the equal-weight average of the air-neighbour values.  Buried interior
        // cells (no air neighbour) are left at the Pass-1 zero state and continue to
        // serve as the Poisson Dirichlet anchor for p_dyn.
        #pragma omp parallel for schedule(static)
        for (int i = 1; i < m.im-1; i++) {
            for (int j = 0; j < m.jm; j++) {
                for (int k = 0; k < m.km; k++) {
                    if (is_land(m.h, i, j, k)) {
                        average_from_air_neighbors(i, j, k);
                    }
                }
            }
        }

        // Pass 3 — copy surface-level (i_mount) values of temperature, microphysics,
        // microphysical sources, and precipitation down to the i=0 reference layer.
        // This keeps the sea-level layer consistent with the lowest atmospheric level
        // so that surface fluxes and output diagnostics see the correct surface state.
        #pragma omp parallel for schedule(static)
        for (int j = 0; j < m.jm; j++) {
            for (int k = 0; k < m.km; k++) {
                int i_mount = m.i_topography[j][k];

                if (is_finite_safe(m.t.x[i_mount][j][k]))
                    m.t.x[0][j][k] = m.t.x[i_mount][j][k];

                m.radiation.x[0][j][k] = m.radiation.x[i_mount][j][k];

                m.p_dyn.x[0][j][k] = m.p_dyn.x[i_mount][j][k];

                m.r_dry.x[0][j][k] = m.r_dry.x[i_mount][j][k];
                m.r_humid.x[0][j][k] = m.r_humid.x[i_mount][j][k];

                m.c.x[0][j][k]     = m.c.x[i_mount][j][k];
                m.cloud.x[0][j][k] = m.cloud.x[i_mount][j][k];
                m.ice.x[0][j][k]   = m.ice.x[i_mount][j][k];
                m.gr.x[0][j][k]    = m.gr.x[i_mount][j][k];

                m.co2.x[0][j][k]   = m.co2.x[i_mount][j][k];

                m.S_c_c.x[0][j][k] = m.S_c_c.x[i_mount][j][k];
                m.S_v.x[0][j][k]   = m.S_v.x[i_mount][j][k];
                m.S_c.x[0][j][k]   = m.S_c.x[i_mount][j][k];
                m.S_i.x[0][j][k]   = m.S_i.x[i_mount][j][k];
                m.S_r.x[0][j][k]   = m.S_r.x[i_mount][j][k];
                m.S_s.x[0][j][k]   = m.S_s.x[i_mount][j][k];

                m.q_v_u.x[0][j][k] = m.q_v_u.x[i_mount][j][k];
                m.q_c_u.x[0][j][k] = m.q_c_u.x[i_mount][j][k];
                m.q_v_d.x[0][j][k] = m.q_v_d.x[i_mount][j][k];

                m.E_u.x[0][j][k] = m.E_u.x[i_mount][j][k];
                m.E_d.x[0][j][k] = m.E_d.x[i_mount][j][k];
                m.D_u.x[0][j][k] = m.D_u.x[i_mount][j][k];
                m.D_d.x[0][j][k] = m.D_d.x[i_mount][j][k];
                m.M_u.x[0][j][k] = m.M_u.x[i_mount][j][k];
                m.M_d.x[0][j][k] = m.M_d.x[i_mount][j][k];

                m.tke.x[0][j][k] = m.tke.x[i_mount][j][k];
                m.dis.x[0][j][k] = m.dis.x[i_mount][j][k];
            }
        }

        // Pass 4 — clip any negative microphysics values to zero across the full
        // domain.  Negative values can arise from the advection scheme near sharp
        // gradients; clamping here prevents unphysical condensate.
        #pragma omp parallel for schedule(static)
        for (int i = 0; i < m.im; i++) {
            for (int j = 0; j < m.jm; j++) {
                for (int k = 0; k < m.km; k++) {
                    if (m.c.x[i][j][k]     < 0.0)  m.c.x[i][j][k]     = 0.0;
                    if (m.cloud.x[i][j][k] < 0.0)  m.cloud.x[i][j][k] = 0.0;
                    if (m.ice.x[i][j][k]   < 0.0)  m.ice.x[i][j][k]   = 0.0;
                    if (m.tke.x[i][j][k]   < 0.0)  m.tke.x[i][j][k]   = 0.0;
                    if (m.dis.x[i][j][k]   < 0.0)  m.dis.x[i][j][k]   = 0.0;
                }
            }
        }

        // Pass 5 — defensive bit-level sanitization of turbulence fields.
        // Under -ffast-math, std::max(0.0, x) and std::clamp(...) are unreliable on
        // NaN inputs (the compiler assumes operands are finite).  TKE reached 7×10²⁰⁰
        // at one polar tropopause cell despite the RK4 cap because a NaN intermediate
        // in prod (computed as cnue * grad·grad − 0.667·tke·grad with grad large at
        // the polar metric singularity) leaked through max(0,...).  Use the same
        // bit-level non-finite check as Paraview_Atm::safe_val to reset non-finite
        // values to a physically sane default, then clip to physical ceilings.
        // This is the last operation in bcSolidGround so all fields are bounded by
        // the time the next RK4 cycle reads them.
        const double tke_max_phys = 1000.0;
        const double tke_max_nd   = tke_max_phys / (m.u_0 * m.u_0);
        const double nue_max      = 1000.0 / (m.u_0 * m.L_atm);
        const double prod_max     = 1.0e4;            // [m²/s³] — generous; real prod ~ 1e0
        const double dis_max      = 1.0e6;            // [m²/s³] — generous for ω*

        // Velocity ceiling — u,v,w had no NaN/clamp guard, yet the 1/sinθ metric
        // drives the radial velocity u into a runaway downdraft at steep high-latitude
        // coasts (~70°W: Baffin 66°N, Antarctic Peninsula 70°S), which then propagates
        // to NaN. 100 m/s is generous (real winds ≲ 110 m/s, model w peaks ~30 m/s) so
        // healthy flow is untouched; this is a stabiliser, not a root-cause cure.
        const double vel_max_phys = 100.0;            // [m/s]
        const double vel_max_nd   = vel_max_phys / m.u_0;

        auto safe_clamp = [](double v, double lo, double hi) -> double {
            std::uint64_t bits;
            std::memcpy(&bits, &v, sizeof(bits));
            if ((bits & 0x7FF0000000000000ULL) == 0x7FF0000000000000ULL) return lo;
            return (v < lo) ? lo : (v > hi) ? hi : v;
        };

        // Symmetric variant for velocities: a non-finite cell relaxes to rest (0),
        // not to the lower bound, and is clipped to ±mag.
        auto safe_clamp_sym = [](double v, double mag) -> double {
            std::uint64_t bits;
            std::memcpy(&bits, &v, sizeof(bits));
            if ((bits & 0x7FF0000000000000ULL) == 0x7FF0000000000000ULL) return 0.0;
            return (v < -mag) ? -mag : (v > mag) ? mag : v;
        };

        #pragma omp parallel for schedule(static)
        for (int i = 0; i < m.im; i++) {
            for (int j = 0; j < m.jm; j++) {
                for (int k = 0; k < m.km; k++) {
                    m.tke.x[i][j][k]        = safe_clamp(m.tke.x[i][j][k],        0.0,        tke_max_nd);
                    m.tken.x[i][j][k]       = safe_clamp(m.tken.x[i][j][k],       0.0,        tke_max_nd);
                    m.dis.x[i][j][k]        = safe_clamp(m.dis.x[i][j][k],        1.0e-10,    dis_max);
                    m.disn.x[i][j][k]       = safe_clamp(m.disn.x[i][j][k],       1.0e-10,    dis_max);
                    m.nue.x[i][j][k]        = safe_clamp(m.nue.x[i][j][k],        0.0,        nue_max);
                    m.prod.x[i][j][k]       = safe_clamp(m.prod.x[i][j][k],       0.0,        prod_max);
                    m.tke_source.x[i][j][k] = safe_clamp(m.tke_source.x[i][j][k], -prod_max,  prod_max);
                    m.dis_source.x[i][j][k] = safe_clamp(m.dis_source.x[i][j][k], -prod_max,  prod_max);
                    m.u.x[i][j][k]          = safe_clamp_sym(m.u.x[i][j][k],          vel_max_nd);
                    m.un.x[i][j][k]         = safe_clamp_sym(m.un.x[i][j][k],         vel_max_nd);
                    m.v.x[i][j][k]          = safe_clamp_sym(m.v.x[i][j][k],          vel_max_nd);
                    m.vn.x[i][j][k]         = safe_clamp_sym(m.vn.x[i][j][k],         vel_max_nd);
                    m.w.x[i][j][k]          = safe_clamp_sym(m.w.x[i][j][k],          vel_max_nd);
                    m.wn.x[i][j][k]         = safe_clamp_sym(m.wn.x[i][j][k],         vel_max_nd);
                }
            }
        }

        cout << "      AGCM: BC_SolidGround ended" << endl;
    }
/*
* 
*/
    // ------------------------------------------------------------------
    void bcRadius()
    {
        // cubic extrapolation: x[a] = x[a+3d] - 3*x[a+2d] + 3*x[a+d]
        // von Neumann:         x[a] = c43*x[a+d] - c13*x[a+2d]

        // Pattern A: ZERO-GRADIENT COPY at both i=0 and i=im-1 (was cubic at both).
        // MC_t/MC_q/MC_v/MC_w are deliberately NOT extrapolated — MoistConvection
        // rhsForcing already writes them (capped) at i=0..im-2, and the 3-point linear
        // extrap x[0]=x[3]-3x[2]+3x[1] overshoots 7× when the capped values oscillate
        // (e.g. +0.01,-0.01,+0.01 → 0.07 at i=0), reintroducing the unbounded MC forcing
        // at the surface and driving the Gulf-of-Alaska coastal velocity runaway.
        Array* both_cubic[] = {
            &m.PrecipitableWaterLocal,
            &m.CoriolisForce, &m.CentrifugalForce, &m.BuoyancyForce, &m.PresGradForce,
            &m.Q_Latent, &m.Q_Sensible,
            &m.S_c_c, &m.S_v, &m.S_c, &m.S_i, &m.S_r, &m.S_s, &m.S_g,
            &m.q_v_u, &m.q_c_u, &m.u_u, &m.v_u, &m.w_u,
            &m.s, &m.s_u, &m.s_d,
            &m.M_u, &m.M_d
        };
        constexpr int n_both = sizeof(both_cubic) / sizeof(both_cubic[0]);

        // Pattern B: von Neumann at i=0, zero-gradient plain copy at i=im-1.
        // The lid used to be the cubic (zero-third-difference) extrapolation
        // x[iml]=x[iml-3]-3x[iml-2]+3x[iml-1]. That preserves the profile's
        // curvature and projects it past the boundary, so a v/w field decaying
        // toward zero at the top overshoots THROUGH zero — the opposite-sign
        // top reflections seen by ~iter 20. The reflection then feeds back via
        // the i=iml-1 radial stencil. A non-amplifying zero-gradient copy
        // (free-slip, no vertical stress at the lid) removes it, mirroring the
        // p_dyn lid fix in PressureSolverAtm.h.
        Array* vn_bot_cubic_top[] = {
            &m.u, &m.v, &m.w, 
            &m.P_conv 
//            &m.tke, &m.dis 
            };
        constexpr int n_vn_cubic = sizeof(vn_bot_cubic_top) / sizeof(vn_bot_cubic_top[0]);

        // Pattern C: cubic at i=0, von Neumann at i=im-1
        // Turbulence scalars use this pattern: at i=im-1 turbulence decays to its
        // background floor, making the profile concave. Cubic extrapolation there
        // amplifies the concavity and produces values much larger than nue[im-2].
        // Von Neumann (zero gradient) is the correct physical BC at the model top.
        // At i=0 cubic is used; apply_wall_bc() reasserts the correct dis value.
//        Array* cubic_bot_vn_top[] = { &m.c, &m.cloud, &m.ice, &m.gr, &m.co2,
        Array* cubic_bot_vn_top[] = {
            &m.cloud, &m.ice, &m.gr, &m.co2,
            &m.c_u, &m.e_d, &m.e_l, &m.e_p, &m.g_p,
            &m.tke, &m.dis, &m.nue, &m.prod, &m.tke_source, &m.dis_source
        };
        constexpr int n_cubic_vn = sizeof(cubic_bot_vn_top) / sizeof(cubic_bot_vn_top[0]);

        const int iml = m.im - 1;

        // Pin t at the lid to its IC value when the snapshot is available;
        // otherwise fall back to a zero-gradient copy (never the old cubic).
        const bool pin_t_top = ((int)m.t_top_init.size() == m.jm);

        #pragma omp parallel for schedule(static)
        for (int j = 0; j < m.jm; j++) {
            for (int k = 0; k < m.km; k++) {

                // t lid (i=im-1): pin to the IC isothermal-floor value. The old
                // cubic extrapolation projected interior curvature onto the lid
                // and amplified it (this stencil's condition number is ~7), so the
                // constant stratospheric top drifted upward within ~20 iters and
                // corrupted the otherwise-steady initial state. Pinning holds it
                // fixed; the upper field then only evolves around orography.
                m.t.x[iml][j][k] = pin_t_top ? m.t_top_init[j][k]
                                             : m.t.x[iml-1][j][k];

                // radiation lid (i=im-1): LEFT ALONE. It used to be set to
                // σ·(t·t_0)⁴ on the reasoning that "radiation.x = σ·(t·t_0)⁴", which
                // stopped being true when MultiLayerRadiation was rewritten as flux
                // sweeps: radiation.x is now the upward long-wave flux at the top of
                // each layer, so radiation.x[im-1] IS the outgoing long-wave flux.
                // Overwriting it with a blackbody at the lid temperature would throw
                // away the one number the column integration exists to produce and
                // put back the artefact that made the OLR a prescribed input.
                //
                // No extrapolation is needed in its place: MultiLayerRadiation writes
                // every level of radiation.x, including the lid, every call.

                // Pattern A — cubic at i=0, ZERO-GRADIENT COPY at the lid (2026-08-19,
                // ported from ATHAD).
                //
                // The lid used to be the same cubic. It is the third time this stencil has had
                // to come off a boundary in this file: v/w "overshoot THROUGH zero", the
                // turbulence scalars "amplify the concavity", p_stat reached -36 hPa in ATHAD
                // (Pattern D below). It was putting NEGATIVE values into fields that cannot be
                // negative — measured here: s_u = -5.01 at the SURFACE, i.e. a parcel at
                // -1020 K. Not one of these fields wants curvature projected past a boundary:
                // convection scalars, microphysical rates, mass fluxes, diagnostic forces. A
                // zero-gradient copy is non-amplifying and cannot manufacture a sign. The lid
                // values are read — rhsForcing's flux divergence at i = im-2 uses
                // M_u/M_d/s/s_u/s_d/q_v_u at i = im-1 — so this is not only cosmetic.
                // AND THE SAME AT i=0, for the same reason and on the same evidence. The
                // surface cubic x[0] = x[3] - 3x[2] + 3x[1] is the identical stencil and it
                // does the identical thing: measured in ATHAD_COND, s_u = -5.012 at 0 m, i.e.
                // an updraft parcel at -1020 K sitting on a 513 K sea. (ATHAD does not show it
                // at the surface, because precompute zeroes s_u there and most of its columns
                // carry no updraft — the same latent defect, exposed by the other fork. A
                // cross-reference is not a check, but a fork IS a second measurement.) This
                // file's own comment already records that this stencil "overshoots 7x" at i=0
                // for the MC_* fields, which is why they are not in this list at all.
                for (int f = 0; f < n_both; f++) {
                    double*** xf = both_cubic[f]->x;
                    xf[0][j][k]   = xf[1][j][k];
                    xf[iml][j][k] = xf[iml-1][j][k];
                }

                // Pattern D — hydrostatic quantities: cubic at i=0, LOG-LINEAR at the lid.
                //
                // p_stat, r_humid and r_dry were on Pattern A (cubic at both ends). The
                // cubic x[iml] = x[iml-3] - 3x[iml-2] + 3x[iml-1] preserves curvature and
                // projects it past the boundary — the same defect this file already
                // documents for v/w ("overshoots THROUGH zero") and for the turbulence
                // scalars ("amplifies the concavity"). Pressure and density are not just
                // concave at the lid, they decay EXPONENTIALLY, so on ATHAD's 300 km shell
                // the cubic drove p_stat to -36 hPa across ~11500 cells at the top level.
                // Harmless on Earth's 16 km shell, where the profile is nearly linear over
                // the last three cells.
                //
                // In log space a hydrostatic profile is straight — exactly straight for an
                // isothermal layer, which is what the top of this column is. So extrapolate
                // the LOGARITHM linearly: x[iml] = x[iml-1]^2 / x[iml-2]. That is the
                // correct continuation of exponential decay and is positive by
                // construction. Falls back to a zero-gradient copy if either sample is
                // non-positive, so it can never introduce the value it exists to prevent.
                {
                    Array* hydro[] = { &m.p_stat, &m.r_humid, &m.r_dry };
                    for (int f = 0; f < 3; f++) {
                        double*** xf = hydro[f]->x;
                        xf[0][j][k] = xf[3][j][k] - 3.0 * xf[2][j][k] + 3.0 * xf[1][j][k];

                        const double x1 = xf[iml-1][j][k];
                        const double x2 = xf[iml-2][j][k];
                        xf[iml][j][k] = (x1 > 0.0 && x2 > 0.0) ? (x1 * x1 / x2) : x1;
                    }
                }

                // Pattern B — von Neumann at i=0, zero-gradient copy at the lid.
                // (u is overwritten with the rigid-lid u=0 just below; P_conv keeps
                // the lid copy; v,w get a taper-to-zero ceiling just after.)
                for (int f = 0; f < n_vn_cubic; f++) {
                    double*** xf = vn_bot_cubic_top[f]->x;
                    xf[0][j][k]   = m.c43 * xf[1][j][k] - m.c13 * xf[2][j][k];
                    xf[iml][j][k] = xf[iml-1][j][k];
                }

                // Horizontal velocity (v,w) taper to a zero-velocity grid ceiling.
                // The free-slip zero-gradient copy above carries the tropopause-jet
                // value straight to i=im-1, so v,w stay finite ("stretched up") at the
                // top; the IC instead ramps v,w to 0 above the tropopause. Restore that
                // by ramping the top three layers smoothly to zero (factor 2/3, 1/3, 0)
                // so the lid is quiet WITHOUT the one-cell shear shock a hard zero would
                // make. This is the horizontal analogue of the u=0 rigid lid below.
                // vn,wn are held in lockstep (as un is for u) so the next RK4 step
                // starts from the tapered ceiling.
                m.v.x[iml][j][k]    = 0.0;          m.w.x[iml][j][k]    = 0.0;
                m.vn.x[iml][j][k]   = 0.0;          m.wn.x[iml][j][k]   = 0.0;
                m.v.x[iml-1][j][k]  *= (1.0/3.0);   m.w.x[iml-1][j][k]  *= (1.0/3.0);
                m.vn.x[iml-1][j][k] *= (1.0/3.0);   m.wn.x[iml-1][j][k] *= (1.0/3.0);
                m.v.x[iml-2][j][k]  *= (2.0/3.0);   m.w.x[iml-2][j][k]  *= (2.0/3.0);
                m.vn.x[iml-2][j][k] *= (2.0/3.0);   m.wn.x[iml-2][j][k] *= (2.0/3.0);

                // Rigid lid on the VERTICAL velocity u at the model top (i=im-1).
                // u is the radial component; air cannot flow through the lid, so the
                // physical top BC is u = 0 — NOT the cubic extrapolation Pattern B just
                // wrote (correct only for the HORIZONTAL v,w, which may have a tropopause
                // jet). The cubic overshoots an increasing u-profile and feeds back via
                // the i=im-2 ∂/∂r stencil, ratcheting upper-level vertical velocity up
                // every step (observed: u at 30–45°N grows 0 → ~17 m/s by iter 300,
                // accelerating, while the zonal jet w is unchanged). The rigid lid removes
                // that feedback AND closes the column mass budget for the all-Neumann
                // pressure Poisson, which otherwise leaves the column-mean vertical
                // velocity an undetermined, drifting constant.
                m.u.x[iml][j][k]  = 0.0;
                m.un.x[iml][j][k] = 0.0;

                // Surface no-penetration BC on the radial velocity u at i=0 (Dirichlet
                // u=0), enforced ALWAYS — not just during the inviscid spin-up.
                // u is the wall-normal component at i=0 (ocean surface above water, solid
                // ground below land i_topography). Tangential v, w keep the Pattern B
                // Neumann (zero-gradient) values — that mirror across the wall is the
                // free-slip tangential BC. un is held in lockstep so the next-time-level
                // field starts from the wall condition.
                //
                // The previous `if (m.inviscid_phase)` guard left u at i=0 unconstrained
                // through the viscous phase. With inviscid_spinup_iters=0 that meant u
                // at the ocean surface was only the Pattern B cubic extrapolation from
                // u[1] and u[2] — no friction, no no-penetration. Coastal pressure
                // gradients then drove u → ±100 m/s at the Gulf-of-Alaska / Tibet / Andes
                // coastal cells (observed iter 240+ → catastrophic t/p/ρ at iter 341,
                // overflow-cascade NaN at iter 359 in v3/v4/v5 runs).
                m.u.x[0][j][k]  = 0.0;
                m.un.x[0][j][k] = 0.0;

                // Pattern C
                for (int f = 0; f < n_cubic_vn; f++) {
                    double*** xf = cubic_bot_vn_top[f]->x;
                    xf[0][j][k]   = xf[3][j][k]     - 3.0 * xf[2][j][k]     + 3.0 * xf[1][j][k];
                    xf[iml][j][k] = m.c43 * xf[iml-1][j][k] - m.c13 * xf[iml-2][j][k];
                }

                // updraft moisture at the surface must be non-negative
                if (m.q_v_u.x[0][j][k] < 0.0) m.q_v_u.x[0][j][k] = 0.0;
                if (m.q_c_u.x[0][j][k] < 0.0) m.q_c_u.x[0][j][k] = 0.0;
            }
        }
    }
/*
* 
*/
    // ------------------------------------------------------------------
    void bcTheta()
    {
        // Pole BC: zero-gradient (copy from first interior cell). At the spherical
        // singularity sin θ → 0, any extrapolation amplifies grid noise:
        //   - (4/3, -1/3) form has amplification factor 4/3
        //   - cubic p[0] = p[3] − 3 p[2] + 3 p[1] has condition number ~7
        // With strong bulk flow these are washed out, but with a weak / spinning-up
        // velocity field the amplified noise dominates and reaches NaN at the pole
        // within ~150 iter (observed at 90°N 0°E, height 0 m). Plain copy gives
        // factor 1.0 — same as an axisymmetric-pole assumption to first order.
        Array* fields_vn[] = {
            &m.t, &m.radiation, &m.p_stat, &m.r_humid, &m.r_dry,
            &m.u, &m.v, &m.w,
            &m.c, &m.cloud, &m.ice, &m.gr,
            &m.PrecipitableWaterLocal, &m.co2,
            &m.tke, &m.dis, &m.nue, &m.prod, &m.tke_source, &m.dis_source
        };
        constexpr int n_vn = sizeof(fields_vn) / sizeof(fields_vn[0]);

        // NAME ONLY, no cubic here: this list gets the same plain copy as fields_vn above
        // (see the loop). Renamed 2026-08-19 from "fields_cubic", which is what it once was —
        // a stale name for a stencil the pole argument above already removed. Checked while
        // repairing the RADIAL cubic (README items 52-53) and found already correct; the name
        // was the only thing left of the defect, and a name is how the next reader decides
        // where to look.
        Array* fields_pole_copy[] = {
            &m.S_c_c, &m.S_v, &m.S_c, &m.S_i, &m.S_r, &m.S_s, &m.S_g,
            &m.q_v_u, &m.q_c_u, &m.u_u, &m.v_u, &m.w_u,
            &m.s, &m.s_u, &m.s_d,
            &m.CoriolisForce, &m.CentrifugalForce, &m.BuoyancyForce, &m.PresGradForce,
            &m.Q_Latent, &m.Q_Sensible,
            &m.c_u, &m.e_d, &m.e_l, &m.e_p, &m.g_p
        };
        constexpr int n_pole_copy = sizeof(fields_pole_copy) / sizeof(fields_pole_copy[0]);

        const int jml = m.jm - 1;

        #pragma omp parallel for schedule(static)
        for (int k = 0; k < m.km; k++) {
            for (int i = 0; i < m.im; i++) {

                for (int f = 0; f < n_vn; f++) {
                    double*** xf = fields_vn[f]->x;
                    xf[i][0][k]   = xf[i][1][k];
                    xf[i][jml][k] = xf[i][jml-1][k];
                }

                for (int f = 0; f < n_pole_copy; f++) {
                    double*** xf = fields_pole_copy[f]->x;
                    xf[i][0][k]   = xf[i][1][k];
                    xf[i][jml][k] = xf[i][jml-1][k];
                }
            }
        }
    }
/*
* 
*/
    // ------------------------------------------------------------------
    void bcPhi()
    {
        // fields: von Neumann extrapolation + periodic averaging
        Array* fields_avg[] = {
            &m.t, &m.radiation, &m.p_stat, &m.r_humid, &m.r_dry,
            &m.u, &m.v, &m.w,
            &m.c, &m.cloud, &m.ice, &m.gr, &m.co2,
            &m.S_c_c, &m.S_v, &m.S_c, &m.S_i, &m.S_r, &m.S_s, &m.S_g,
            &m.q_v_u, &m.q_c_u, &m.u_u, &m.v_u, &m.w_u,
            &m.s, &m.s_u, &m.s_d,
            &m.CoriolisForce, &m.CentrifugalForce, &m.PresGradForce, &m.BuoyancyForce,
            &m.Q_Latent, &m.Q_Sensible,
            &m.c_u, &m.e_d, &m.e_l, &m.e_p, &m.g_p
//            &m.tke, &m.dis, &m.nue, &m.prod, &m.tke_source, &m.dis_source
        };
        constexpr int n_avg = sizeof(fields_avg) / sizeof(fields_avg[0]);

        // fields: von Neumann only (no periodic averaging)
        Array* fields_extrap[] = { 
            &m.PrecipitableWaterLocal, 
            &m.tke, &m.dis, &m.nue, &m.prod, &m.tke_source, &m.dis_source
        };
        constexpr int n_extrap = sizeof(fields_extrap) / sizeof(fields_extrap[0]);

        #pragma omp parallel for schedule(static)
        for (int i = 0; i < m.im; i++) {
            for (int j = 0; j < m.jm; j++) {

                for (int f = 0; f < n_avg; f++) {
                    double** xij = fields_avg[f]->x[i];
                    double v0   = m.c43 * xij[j][1]      - m.c13 * xij[j][2];
                    double vend = m.c43 * xij[j][m.km-2] - m.c13 * xij[j][m.km-3];
                    xij[j][0] = xij[j][m.km-1] = (v0 + vend) * 0.5;
                }

                for (int f = 0; f < n_extrap; f++) {
                    double** xij = fields_extrap[f]->x[i];
                    xij[j][0]      = m.c43 * xij[j][1]      - m.c13 * xij[j][2];
                    xij[j][m.km-1] = m.c43 * xij[j][m.km-2] - m.c13 * xij[j][m.km-3];
                }
            }
        }

        // ------------------------------------------------------------------
        // Seam zonal (φ) damping.
        //
        // k=0 and k=km-1 are the SAME physical longitude (0°≡360°) and are NOT
        // evolved by the RK4 (its φ-loop runs k=1..km-2); the reconstruction above
        // pins them to 0.5·(x[1]+x[km-2]).  Because that slaves the seam value to its
        // own neighbours, the discrete d²/dφ² self-damping at the seam-adjacent cells
        // k=1 and k=km-2 drops from −2 to −1.5 — a 25% loss of numerical zonal
        // diffusion exactly at the seam.  Combined with the 1/sin²θ metric and a
        // coastline this leaves an under-damped zonal mode that runs the velocity
        // away: observed at 50°N/1°E (k=1), the ocean cell next to the European coast
        // (land at 2°E) and adjacent to the Greenwich seam.  Its exponential growth is
        // then mirrored into k=0 by the average (k=0 ≈ 0.5·u[1]), which is why the
        // blow-up *looks* like a k=0 problem.  The polar zonal filter does 0 passes at
        // 50°N (sin_ref/sinθ ratio rounds to 0 there), so it cannot catch this.
        //
        // One explicit 1-2-1 Shapiro pass (coeff 0.25 fully removes the 2Δφ mode in a
        // single step) across the three seam cells {km-2, 0≡km-1, 1} restores the lost
        // damping locally, leaving the rest of the field untouched.  Periodic
        // neighbours are used with no-flux at solid cells (a land neighbour relaxes the
        // cell toward the seam instead of dragging it to zero).  Both the time-level-n
        // state (un,vn,wn — what the next RK4 integrates from; storeIntermediateData3D
        // copies u→un, so these hold the previous result) and the current field
        // (u,v,w — what the RHS reads for φ-derivatives) are damped.
        {
            const double seam_coeff = 0.25;
            const int km2 = m.km - 2;
            const int km3 = m.km - 3;

            auto smooth_seam = [&](Array& f, int i, int j) {
                // Snapshot the three seam cells (k=0 and k=km-1 are one point).
                const double f_km2 = f.x[i][j][km2];
                const double f_0   = f.x[i][j][0];
                const double f_1   = f.x[i][j][1];

                const bool air_km2 = is_air(m.h, i, j, km2);
                const bool air_0   = is_air(m.h, i, j, 0);
                const bool air_1   = is_air(m.h, i, j, 1);

                // Periodic neighbours, no-flux (substitute the cell value) at solid.
                const double w_km2 = is_air(m.h, i, j, km3) ? f.x[i][j][km3] : f_km2;
                const double e_km2 = air_0 ? f_0   : f_km2;   // east of km-2 is the seam
                const double w_0   = air_km2 ? f_km2 : f_0;   // west of seam is km-2
                const double e_0   = air_1 ? f_1   : f_0;     // east of seam is k=1
                const double w_1   = air_0 ? f_0   : f_1;     // west of k=1 is the seam
                const double e_1   = is_air(m.h, i, j, 2) ? f.x[i][j][2] : f_1;

                if (air_km2)
                    f.x[i][j][km2] = f_km2 + seam_coeff * (w_km2 - 2.0 * f_km2 + e_km2);
                if (air_0)
                    f.x[i][j][0] = f.x[i][j][m.km-1] =
                        f_0 + seam_coeff * (w_0 - 2.0 * f_0 + e_0);
                if (air_1)
                    f.x[i][j][1] = f_1 + seam_coeff * (w_1 - 2.0 * f_1 + e_1);
            };

            // 2 passes — one 1-2-1 fully damps a 2Δφ mode in steady state, but the
            // upper-troposphere storm-track jet crossing the seam at 62°N (i~34, the
            // iter-349 dry-NaN seed) regenerates seam energy each step faster than a
            // single pass can absorb. A second pass extends the damping reach to ~4Δ.
            #pragma omp parallel for schedule(static)
            for (int i = 0; i < m.im; i++) {
                for (int j = 0; j < m.jm; j++) {
                    for (int p = 0; p < 2; p++) {
                        smooth_seam(m.u, i, j);   smooth_seam(m.un, i, j);
                        smooth_seam(m.v, i, j);   smooth_seam(m.vn, i, j);
                        smooth_seam(m.w, i, j);   smooth_seam(m.wn, i, j);
                    }
                }
            }
        }
    }
/*
*
*/
    // ------------------------------------------------------------------
    void bcVelSurfSur()
    {
        using namespace std;
        cout << endl << "      AGCM: BC_vel_surf_sur" << endl;


        constexpr double coeff  = 0.01;                                 // represents the flow better
        constexpr double coeff5 = 0.9;                                  // represents the flow better

        #pragma omp parallel for schedule(static)
        for (int i = 1; i < m.im-1; i++) {
            for (int j = 1; j < m.jm-1; j++) {
                for (int k = 1; k < m.km-1; k++) {

                    bool land_ijk = is_land(m.h, i, j, k);

                    // i-direction (radial: land→air outward)
                    if (i < m.im-2 && land_ijk && is_air(m.h, i+1, j, k)) {
                        m.u.x[i][j][k]   = 0.0;
                        m.u.x[i+1][j][k] = 0.0;

                        m.v.x[i][j][k]    = 0.0;
                        m.v.x[i+1][j][k] *= coeff;
                        m.v.x[i+2][j][k] *= coeff5;

                        m.w.x[i][j][k]    = 0.0;
                        m.w.x[i+1][j][k] *= coeff;
                        m.w.x[i+2][j][k] *= coeff5;
                    }

                    // j-direction (meridional: land→air south/north)
                    if (land_ijk) {
                        if (j < m.jm-2 && is_air(m.h, i, j+1, k) && is_air(m.h, i, j+2, k)) {  // south
                            m.u.x[i][j][k]    = 0.0;
                            m.u.x[i][j+1][k] *= coeff;
                            m.u.x[i][j+2][k] *= coeff5;

                            m.v.x[i][j][k]    = 0.0;
                            m.v.x[i][j+1][k]  = 0.0;

                            m.w.x[i][j][k]    = 0.0;
                            m.w.x[i][j+1][k] *= coeff;
                            m.w.x[i][j+2][k] *= coeff5;
                        }

                        if (j >= 2 && is_air(m.h, i, j-1, k) && is_air(m.h, i, j-2, k)) {  // north
                            m.u.x[i][j][k]    = 0.0;
                            m.u.x[i][j-1][k] *= coeff;
                            m.u.x[i][j-2][k] *= coeff5;

                            m.v.x[i][j][k]    = 0.0;
                            m.v.x[i][j-1][k]  = 0.0;

                            m.w.x[i][j][k]    = 0.0;
                            m.w.x[i][j-1][k] *= coeff;
                            m.w.x[i][j-2][k] *= coeff5;
                        }
                    }

                    // k-direction (zonal: land→air east/west)
                    if (land_ijk) {
                        if (k < m.km-2 && is_air(m.h, i, j, k+1) && is_air(m.h, i, j, k+2)) {  // east
                            m.u.x[i][j][k]    = 0.0;
                            m.u.x[i][j][k+1] *= coeff;
                            m.u.x[i][j][k+2] *= coeff5;

                            m.v.x[i][j][k]    = 0.0;
                            m.v.x[i][j][k+1] *= coeff;
                            m.v.x[i][j][k+2] *= coeff5;

                            m.w.x[i][j][k]    = 0.0;
                            m.w.x[i][j][k+1]  = 0.0;
                        }

                        if (k >= 2 && is_air(m.h, i, j, k-1) && is_air(m.h, i, j, k-2)) {  // west
                            m.u.x[i][j][k]    = 0.0;
                            m.u.x[i][j][k-1] *= coeff;
                            m.u.x[i][j][k-2] *= coeff5;

                            m.v.x[i][j][k]    = 0.0;
                            m.v.x[i][j][k-1] *= coeff;
                            m.v.x[i][j][k-2] *= coeff5;

                            m.w.x[i][j][k]    = 0.0;
                            m.w.x[i][j][k-1]  = 0.0;
                        }
                    }
                }
            }
        }
        cout << "      AGCM: BC_vel_surf_sur ended" << endl;
    }
/*
*
*/
    // ------------------------------------------------------------------
    // Neumann (zero normal gradient) boundary condition for all scalar
    // 3-D fields at land-air interfaces.  For each interior land cell the
    // value is extrapolated from the two nearest air cells in the normal
    // direction: x[land] = c43*x[n1] - c13*x[n2].  Velocity components
    // from MoistConvection (u_u, v_u, w_u, u_d, v_d, w_d, MC_v, MC_w)
    // are excluded; the main-model velocities are handled by bcVelSurfSur.
    // At corners where multiple directions border air, the vertical (i)
    // direction is preferred.
    void bcScalarSurfSur()
    {
        using namespace std;
        cout << endl << "      AGCM: BC_scalar_surf_sur" << endl;

        // S-functions (S_v, S_c, S_i, S_r, S_s, S_g, S_c_c), precipitation
        // (P_rain, P_snow, Precipitation), and p_dyn are excluded: bcSolidGround
        // zeros them at land cells and projects x[i_mount] to i=0, so Neumann
        // extrapolation here would only produce spuriously high land-cell values.
        Array* scalars[] = {
            &m.t, &m.radiation,
            &m.c, &m.cloud, &m.ice, &m.gr, &m.co2,
            &m.p_stat,
            &m.r_dry, &m.r_humid,
            &m.Q_Latent, &m.Q_Sensible,
            &m.BuoyancyForce, &m.CoriolisForce, &m.CentrifugalForce, &m.PresGradForce,
            &m.PrecipitableWaterLocal,
            &m.M_u,   &m.M_d,
            &m.MC_t,  &m.MC_q,
            &m.s,     &m.s_u,   &m.s_d,
            &m.q_v_u, &m.q_v_d, &m.q_c_u,
            &m.c_u,   &m.e_d,   &m.e_l,  &m.e_p,  &m.g_p,
            &m.E_u,   &m.D_u,   &m.E_d,  &m.D_d
        };
        constexpr int ns = sizeof(scalars) / sizeof(scalars[0]);

        #pragma omp parallel for schedule(static)
        for (int i = 1; i < m.im-1; i++) {
            for (int j = 1; j < m.jm-1; j++) {
                for (int k = 1; k < m.km-1; k++) {

                    if (!is_land(m.h, i, j, k)) continue;

                    // ---- i-direction (vertical, preferred at corners) ----
                    if (is_air(m.h, i+1, j, k)) {
                        if (i + 2 < m.im) {
                            for (int f = 0; f < ns; f++) {
                                double*** x = scalars[f]->x;
                                x[i][j][k] = m.c43 * x[i+1][j][k] - m.c13 * x[i+2][j][k];
                            }
                        } else {
                            for (int f = 0; f < ns; f++)
                                scalars[f]->x[i][j][k] = scalars[f]->x[i+1][j][k];
                        }
                        continue;                                    // i preferred: skip j and k
                    }

                    // ---- j-direction ----
                    if (j < m.jm-2 && is_air(m.h, i, j+1, k) && is_air(m.h, i, j+2, k)) {
                        for (int f = 0; f < ns; f++) {
                            double*** x = scalars[f]->x;
                            x[i][j][k] = m.c43 * x[i][j+1][k] - m.c13 * x[i][j+2][k];
                        }
                    } else if (j >= 2 && is_air(m.h, i, j-1, k) && is_air(m.h, i, j-2, k)) {
                        for (int f = 0; f < ns; f++) {
                            double*** x = scalars[f]->x;
                            x[i][j][k] = m.c43 * x[i][j-1][k] - m.c13 * x[i][j-2][k];
                        }
                    }

                    // ---- k-direction ----
                    if (k < m.km-2 && is_air(m.h, i, j, k+1) && is_air(m.h, i, j, k+2)) {
                        for (int f = 0; f < ns; f++) {
                            double*** x = scalars[f]->x;
                            x[i][j][k] = m.c43 * x[i][j][k+1] - m.c13 * x[i][j][k+2];
                        }
                    } else if (k >= 2 && is_air(m.h, i, j, k-1) && is_air(m.h, i, j, k-2)) {
                        for (int f = 0; f < ns; f++) {
                            double*** x = scalars[f]->x;
                            x[i][j][k] = m.c43 * x[i][j][k-1] - m.c13 * x[i][j][k-2];
                        }
                    }
                }
            }
        }

        #pragma omp parallel for schedule(static)
        for (int j = 0; j < m.jm; j++) {
            for (int k = 0; k < m.km; k++) {
                int i_mount = m.i_topography[j][k];

                if (is_finite_safe(m.t.x[i_mount][j][k]))
                    m.t.x[0][j][k] = m.t.x[i_mount][j][k];
                m.radiation.x[0][j][k] = m.radiation.x[i_mount][j][k];
                m.c.x[0][j][k]     = m.c.x[i_mount][j][k];
                m.cloud.x[0][j][k] = m.cloud.x[i_mount][j][k];
                m.ice.x[0][j][k]   = m.ice.x[i_mount][j][k];
                m.gr.x[0][j][k]    = m.gr.x[i_mount][j][k];

                m.S_c_c.x[0][j][k] = m.S_c_c.x[i_mount][j][k];
                m.S_v.x[0][j][k]   = m.S_v.x[i_mount][j][k];
                m.S_c.x[0][j][k]   = m.S_c.x[i_mount][j][k];
                m.S_i.x[0][j][k]   = m.S_i.x[i_mount][j][k];
                m.S_r.x[0][j][k]   = m.S_r.x[i_mount][j][k];
                m.S_s.x[0][j][k]   = m.S_s.x[i_mount][j][k];

                m.P_rain.x[0][j][k] = m.P_rain.x[i_mount][j][k];
                m.P_snow.x[0][j][k] = m.P_snow.x[i_mount][j][k];

                m.p_dyn.x[0][j][k] = m.p_dyn.x[i_mount][j][k];
            }
        }

        cout << "      AGCM: BC_scalar_surf_sur ended" << endl;
    }
/*
* 
*/
    // ------------------------------------------------------------------
    void coastalCurrents()
    {
        using namespace std;
        cout << endl << "      AGCM: CoastalCurrents" << endl;

        constexpr int l_end = 20;

        #pragma omp parallel for schedule(static)
        for (int i = 0; i < m.im-1; i++) {
            for (int j = 1; j < m.jm-1; j++) {
                for (int k = 1; k < m.km-1; k++) {

                    if (!is_land(m.h, i, j, k)) continue;

                    if (is_land(m.h, i, j, k-1) && is_air(m.h, i, j, k+1)) {  // east coast
                        int l_max = std::min(l_end, m.km - 1 - k);
                        for (int l = 1; l <= l_max; l++) {
                            m.v.x[i][j][k+l] = -m.v.x[i][j][k+l];
                            m.w.x[i][j][k+l] = -m.w.x[i][j][k+l];
                        }
                    }

                    if (is_land(m.h, i, j, k+1) && is_air(m.h, i, j, k-1)) {  // west coast
                        int l_max = std::min(l_end, k);
                        for (int l = 1; l <= l_max; l++) {
                            m.v.x[i][j][k-l] = -m.v.x[i][j][k-l];
                            m.w.x[i][j][k-l] = -m.w.x[i][j][k-l];
                        }
                    }
                }
            }
        }
        cout << "      AGCM: CoastalCurrents ended" << endl;
    }

private:
    cAtmosphereModel& m;
};
