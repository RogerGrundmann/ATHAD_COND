#!/usr/bin/env python3
"""Does the prescribed multi-cell structure SURVIVE, or collapse to one cell?

Counts sign bands of Psi across the northern hemisphere, by level and iteration.
5 bands = the prescribed n_cells_hemisphere is intact; 1 = collapsed to a single
overturning. Reported against |Psi| so bands carrying no mass are not counted.

usage: survival.py <outdir> [<outdir> ...]
"""
import re, sys, os

# ATHAD_COND grid: im = 61 (ATHAD is 41). ATM_IM overrides for an off-default grid.
IM = int(os.environ.get("ATM_IM", 61))

def load(fn, im=IM, jm=181):
    txt = open(fn).read().split('\n'); f = {}; i = 0
    while i < len(txt):
        m = re.match(r'SCALARS (\S+) float', txt[i].strip())
        if m:
            n = m.group(1); i += 2
            v = [float(x) for x in txt[i:i+im*jm]]
            f[n] = [[v[a*jm+b] for b in range(jm)] for a in range(im)]
            i += im*jm; continue
        i += 1
    return f

def bands(P, i, mx, thr=1e-2):
    s = [(1 if P[i][j] > 0 else -1) for j in range(1, 90) if abs(P[i][j]) > thr*mx]
    if not s: return 0
    return 1 + sum(1 for a, b in zip(s, s[1:]) if a != b)

for d in sys.argv[1:]:
    print(f"\n=== {d} ===")
    its = sorted(int(m.group(1)) for m in
                 (re.search(r'_(\d+)\.vtk$', f) for f in os.listdir(d) if 'zonal' in f) if m)
    hdr = None
    for it in its:
        fn = [f for f in os.listdir(d) if 'zonal' in f and f.endswith(f'_{it}.vtk')][0]
        f = load(os.path.join(d, fn)); P = f['PsiMerid']; H = f['height']
        mx = max(abs(v) for r in P for v in r)
        if mx == 0:
            print(f"  iter {it:3d}: Psi identically 0"); continue
        levs = [i for i in range(0, IM - 4, max(1, (IM - 4) // 9))]
        if hdr is None:
            hdr = "  iter  " + ' '.join(f'{H[i][90]:>6.0f}km' for i in levs)
            print(hdr)
        print(f"  {it:4d}  " + ' '.join(f'{bands(P,i,mx):>8d}' for i in levs))
