import re,sys
def last(path):
    t=open(path,errors='replace').read()
    d={}
    for m in re.finditer(r'(max|min) ([A-Za-z_0-9][A-Za-z_0-9 \-/]*?)\s*\.* *=\s*(-?[\d.]+(?:e[-+]?\d+)?)',t):
        d[(m.group(2).strip(),m.group(1))]=float(m.group(3))   # last occurrence wins
    return d
a=last(sys.argv[1]+'/run.log'); b=last(sys.argv[2]+'/run.log')
keys=[k for k in a if k in b]; same=[k for k in keys if a[k]==b[k]]
thr=float(sys.argv[3]) if len(sys.argv)>3 else 0.0
print(f"{sys.argv[1]} vs {sys.argv[2]}: {len(same)} of {len(keys)} printed extrema identical at the last print")
for k in keys:
    if a[k]!=b[k]:
        r=(b[k]-a[k])/abs(a[k])*100 if a[k] else float('inf')
        if abs(r)>=thr: print(f"   {k[1]} {k[0]:28s} {a[k]:18.6f} -> {b[k]:18.6f}  {r:+9.3f} %")
