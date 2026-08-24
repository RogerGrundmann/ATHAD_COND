#!/usr/bin/env python3
"""Psi closure metric: Psi(ground) must be 0 for a mass-conserving zonal-mean flow.

Both boundaries are closed (u == 0 at i=0 measured; Psi == 0 at the lid by construction),
so the column-integrated meridional mass flux across every latitude circle must vanish.
Psi(ground)/Psi(interior) is how badly it doesn't, in units of the circulation we want to see.

The interior is selected by PHYSICAL HEIGHT, not level index: a zeta scan moves the heights
under a fixed im, so a fixed index would compare different altitudes between arms.

usage: psicheck.py <label>=<vtk> [...]
"""
import re, sys, os

# ATHAD_COND grid: im = 61 (ATHAD is 41). ATM_IM overrides for an off-default grid.
IM = int(os.environ.get("ATM_IM", 61))
JM = 181

def load(fn, im=IM, jm=JM):
    txt = open(fn).read().split('\n')
    fields, i = {}, 0
    while i < len(txt):
        m = re.match(r'SCALARS (\S+) float', txt[i].strip())
        if m:
            name = m.group(1); i += 2
            vals = [float(x) for x in txt[i:i+im*jm]]
            fields[name] = [[vals[a*jm+b] for b in range(jm)] for a in range(im)]
            i += im*jm; continue
        i += 1
    return fields

H_MIN = 20.0   # km: floor of the "interior" circulation

def metrics(fn):
    f = load(fn); P = f['PsiMerid']; U = f['u-Component']; H = f['height']  # height in km
    lo = [i for i in range(IM) if H[i][90] >= H_MIN]
    ground = max(abs(P[0][j]) for j in range(JM))
    inter  = max(abs(P[i][j]) for i in lo for j in range(JM))
    usurf  = max(abs(U[0][j]) for j in range(JM))
    return ground, inter, usurf, lo[0], H[lo[0]][90]

print(f"{'arm':>14} {'Psi(ground)':>12} {'Psi(>20km)':>12} {'ratio':>7} {'|u|surf':>8} {'i>=20km':>10}")
for arg in sys.argv[1:]:
    lab, fn = arg.split('=', 1)
    g, a, u, i0, h0 = metrics(fn)
    print(f"{lab:>14} {g:12.4e} {a:12.4e} {g/a:7.3f} {u:8.1e}   {i0:2d}@{h0:5.1f}km")
