# field-level comparison of two runs' section VTKs (zonal, longal, radial) at every checkpoint both have:
#   python3 cmp_fields.py <run a> <run b> [field ...]
# per field and checkpoint: largest |b-a| over the three sections, its rms, both relative to the range of the field in run a
import sys,glob,os,re
import numpy as np
FIELDS=['u-Component','v-Component','w-Component','Temperature','WaterVapour','CloudWater','CloudIce',
        'Precipitation','PressureDynamic','Radiation']
def read(path):
    L=open(path).read().split('\n'); d={}; n=0; i=0
    while i<len(L):
        s=L[i]
        if s.startswith('POINT_DATA'): n=int(s.split()[1])
        if s.startswith('SCALARS') and n:
            d[s.split()[1]]=np.array(L[i+2:i+2+n],dtype=float); i+=n+2; continue
        i+=1
    return d
a,b=sys.argv[1],sys.argv[2]; fields=sys.argv[3:] or FIELDS
its=sorted({int(re.search(r'_(\d+)\.vtk$',f).group(1)) for f in glob.glob(a+'/out/*_zonal_*_*.vtk')})
res={f:[] for f in fields}
for it in its:
    acc={f:[0.0,0.0,0,0.0] for f in fields}   # max|d|, sum d^2, count, range
    for fa in glob.glob(a+'/out/*_%d.vtk'%it):
        fb=b+'/out/'+os.path.basename(fa)
        if not os.path.exists(fb): continue
        A=read(fa); Bd=read(fb)
        for f in fields:
            if f not in A or f not in Bd: continue
            d=Bd[f]-A[f]; x=acc[f]
            x[0]=max(x[0],np.abs(d).max()); x[1]+=float((d*d).sum()); x[2]+=d.size; x[3]=max(x[3],A[f].max()-A[f].min())
    for f in fields:
        x=acc[f]; res[f].append((x[0],(x[1]/x[2])**0.5 if x[2] else 0.0,x[3]))
print('largest |difference| as a fraction of the field range, per checkpoint:')
print('%-16s'%'iteration'+''.join('%10d'%i for i in its))
for f in fields: print('%-16s'%f+''.join('%10.1e'%(m/r if r else 0.0) for m,_,r in res[f]))
print('at the last checkpoint: field, range, largest |difference|, rms difference')
for f in fields:
    m,s,r=res[f][-1]; print('  %-16s %14.6g %14.6g %14.6g'%(f,r,m,s))
