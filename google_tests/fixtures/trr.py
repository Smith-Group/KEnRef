import struct, numpy as np
def read_trr(path):
    d=open(path,'rb').read(); o=0
    def i4():
        nonlocal o; v=struct.unpack_from('>i',d,o)[0]; o+=4; return v
    magic=i4(); assert magic==1993, f"not a TRR (magic {magic})"
    slen=i4(); 
    # gromacs writes string as: int len, int len, chars padded  (XDR string of the version)
    o2=o
    nchars=struct.unpack_from('>i',d,o)[0]; o+=4
    o+=(nchars+3)//4*4
    sizes={}
    for k in ('ir','e','box','vir','pres','top','sym','x','v','f'):
        sizes[k]=i4()
    natoms=i4(); step=i4(); nre=i4()
    isdouble = sizes['box']==9*8 if sizes['box'] else sizes['x']==natoms*3*8
    rt='>d' if isdouble else '>f'; rs=8 if isdouble else 4
    t=struct.unpack_from(rt,d,o)[0]; o+=rs
    lam=struct.unpack_from(rt,d,o)[0]; o+=rs
    box=None
    if sizes['box']:
        box=np.array(struct.unpack_from('>'+('d' if isdouble else 'f')*9,d,o)).reshape(3,3); o+=sizes['box']
    for k in ('vir','pres'):
        o+=sizes[k]
    def arr(n):
        nonlocal o
        if not sizes[n]: return None
        a=np.array(struct.unpack_from('>'+('d' if isdouble else 'f')*(natoms*3),d,o)).reshape(natoms,3)
        o+=sizes[n]; return a
    x=arr('x'); v=arr('v'); f=arr('f')
    return dict(natoms=natoms,step=step,t=t,double=isdouble,box=box,x=x,v=v,f=f)
if __name__=='__main__':
    import sys
    r=read_trr(sys.argv[1])
    print(f"natoms={r['natoms']}  step={r['step']}  t={r['t']}  double={r['double']}")
    print("box (nm):\n", r['box'])
    print("box diagonal (nm), repr:", [repr(float(r['box'][i][i])) for i in range(3)])
    print("first 3 atoms (nm):\n", r['x'][:3])
# Minimal TRR reader, enough to pull one frame's box and positions without a GROMACS install.
# TRR is XDR (big-endian); precision is inferred from the recorded box/x array sizes.
