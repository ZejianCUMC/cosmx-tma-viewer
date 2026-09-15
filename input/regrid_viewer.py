#!/usr/bin/env python
# Re-arrange TMA cores of a sample viewer into a tidy N-column grid (default 3),
# discarding the arbitrary TMA core positions. Moves each core's cells AND its mosaic
# tile to the new grid slot. Bakes the render y-flip so cells+mosaic share one space.
import json, re, sys, base64, argparse
from io import BytesIO
import numpy as np, pandas as pd
from scipy.spatial import cKDTree
try:
    from PIL import Image
    Image.MAX_IMAGE_PIXELS=None
except Exception:
    Image=None

def extract_data(txt):
    m=re.search(r"const DATA\s*=\s*", txt); s=m.end(); d=0; i=s
    while i<len(txt):
        c=txt[i]
        if c=='{': d+=1
        elif c=='}':
            d-=1
            if d==0: i+=1; break
        i+=1
    return s,i,json.loads(txt[s:i])

def cores_from_fovs(x,y,fov):
    # Group a sample's FOVs into TMA cores. CosMx numbers each core's FOVs CONSECUTIVELY in
    # raster order (3 cols), so consecutive PRESENT fov numbers are always spatially adjacent
    # WITHIN a core. A core can contain a tissue-free FOV (no cells -> absent from the data);
    # that leaves a ~2x-pitch spatial hole that a pure centroid-distance union splits in two
    # (the historical bug: <unit> fov 124 empty split {122,123} from {125,126,127,128}).
    # Fix = union in two passes:
    #   (1) spatial adjacency for ALL pairs at the original 1.8x-pitch reach (unchanged) -> the
    #       new partition can only COARSEN the old one, never introduce a new split; plus
    #   (2) number-proximity bridges: join near-in-number FOVs (raster-consecutive, <=2 empty FOVs
    #       apart) up to 2.75x pitch, spanning an empty-FOV hole / row-wrap while staying BELOW the
    #       ~2.85-3.1x-pitch spacing that separates distinct TMA cores (so cores are NOT merged).
    #       2.75x is centred between the largest true intra-core row-wrap seen (<unit> 215-216, 2.68x)
    #       and the tightest true inter-core gap (<unit> 330-331, 2.82x); ~25um margin either side.
    cen=pd.DataFrame({"x":x,"y":y,"fov":fov}).groupby("fov")[["x","y"]].mean()
    C=cen[["x","y"]].values; fids=list(cen.index)
    if len(C)==1: return {fids[0]:0}, cen
    tree=cKDTree(C); dnn,_=tree.query(C,k=2); pitch=np.median(dnn[:,1])
    near=max(pitch*1.8, pitch+200)      # spatial adjacency, all pairs (original threshold)
    bridge=max(pitch*2.75, pitch+350)   # extra reach for near-in-number fov pairs only
    parent=list(range(len(C)))
    def find(a):
        while parent[a]!=a: parent[a]=parent[parent[a]]; a=parent[a]
        return a
    def union(a,b): parent[find(a)]=find(b)
    for a,b in tree.query_pairs(near): union(a,b)                    # (1) dense core interior
    fs=np.asarray(fids)                                              # fov numbers, ascending
    for a in range(len(fs)):                                         # (2) raster number-proximity bridges
        for b in range(a+1,len(fs)):
            if fs[b]-fs[a] > 3: break                                # only <=2 tissue-free FOVs between them
            if float(np.hypot(*(C[a]-C[b]))) <= bridge: union(a,b)   # ...and within the intra-core reach
    comp={}
    for i in range(len(C)):
        r=find(i); comp.setdefault(r,[]).append(i)
    core_of_fov={};
    for ci,(r,members) in enumerate(comp.items()):
        for mi in members: core_of_fov[fids[mi]]=ci
    return core_of_fov, cen

def crop(img, box):  # box world (x0,y0,x1,y1) -> PIL crop using ppw
    pass

def regrid(path, out=None, ncol=3, pad=0.06, qc=None):
    txt=open(path,encoding="utf-8").read(); s,e,D=extract_data(txt)
    x=np.array(D["x"],float); y=np.array(D["y"],float); fov=np.array(D.get("fov",[0]*len(x)))
    x0,y0=x.copy(),y.copy()
    if D.get("global"):
        a,b=y.min(),y.max(); y=b-(y-a); D["global"]=False   # bake flip -> display space
    core_of_fov,cen=cores_from_fovs(x,y,fov)
    core=np.array([core_of_fov[f] for f in fov])
    ncore=len(set(core))
    # per-core bbox (padded) from cells
    info={}
    for c in range(ncore):
        m=core==c; cx0,cx1=x[m].min(),x[m].max(); cy0,cy1=y[m].min(),y[m].max()
        w=cx1-cx0; h=cy1-cy0; px=w*pad; py=h*pad
        info[c]=dict(x0=cx0-px,x1=cx1+px,y0=cy0-py,y1=cy1+py,cx=(cx0+cx1)/2,cy=(cy0+cy1)/2,w=w*(1+2*pad),h=h*(1+2*pad))
    cellW=max(v["w"] for v in info.values()); cellH=max(v["h"] for v in info.values())
    # order cores in reading order (top rows first, then left->right)
    order=sorted(info, key=lambda c:(round(info[c]["cy"]/(cellH*0.6)), info[c]["cx"]))
    nrow=int(np.ceil(ncore/ncol))
    # target slot center for each core
    slot={}
    for k,c in enumerate(order):
        col=k%ncol; row=k//ncol
        slot[c]=((col+0.5)*cellW,(row+0.5)*cellH)
    # move cells
    nx=x.copy(); ny=y.copy()
    for c in range(ncore):
        m=core==c; tx,ty=slot[c]; nx[m]+= tx-info[c]["cx"]; ny[m]+= ty-info[c]["cy"]
    newext=[0.0,0.0,ncol*cellW,nrow*cellH]
    # re-tile mosaic
    if Image is not None and D.get("rawImage"):
        e0,e1,e2,e3=D["rawExtent"]; im=Image.open(BytesIO(base64.b64decode(D["rawImage"].split(",",1)[1]))).convert("RGB")
        IW,IH=im.size; ppx=IW/(e2-e0); ppy=IH/(e3-e1)
        canvas=Image.new("RGB",(max(1,int(round(newext[2]*ppx))),max(1,int(round(newext[3]*ppy)))),(0,0,0))
        for c in range(ncore):
            v=info[c]
            sx0=int(round((v["x0"]-e0)*ppx)); sx1=int(round((v["x1"]-e0)*ppx))
            sy0=int(round((v["y0"]-e1)*ppy)); sy1=int(round((v["y1"]-e1)*ppy))
            sx0=max(0,sx0); sy0=max(0,sy0); sx1=min(IW,sx1); sy1=min(IH,sy1)
            if sx1<=sx0 or sy1<=sy0: continue
            tile=im.crop((sx0,sy0,sx1,sy1))
            tx,ty=slot[c]; nwx0=(tx-v["w"]/2); nwy0=(ty-v["h"]/2)
            dx=int(round(nwx0*ppx)); dy=int(round(nwy0*ppy))
            canvas.paste(tile,(dx,dy))
        buf=BytesIO(); canvas.save(buf,format="JPEG",quality=84)
        D["rawImage"]="data:image/jpeg;base64,"+base64.b64encode(buf.getvalue()).decode()
    D["rawExtent"]=newext; D["x"]=[round(float(v),1) for v in nx]; D["y"]=[round(float(v),1) for v in ny]
    newjson=json.dumps(D,separators=(",",":")); out=out or path
    open(out,"w",encoding="utf-8").write(txt[:s]+newjson+txt[e:])
    print(f"[regrid] {path.split('/')[-1]}: {ncore} cores -> {nrow}x{ncol} grid  cell {cellW:.0f}x{cellH:.0f}um")
    if qc:
        import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
        fig,ax=plt.subplots(1,2,figsize=(13,5))
        yd=y0.max()+y0.min()-y0
        ax[0].scatter(x0,yd,s=1,c=core,cmap="tab10"); ax[0].set_title(f"BEFORE ({ncore} cores, TMA positions)"); ax[0].set_aspect("equal")
        ax[1].scatter(nx,ny,s=1,c=core,cmap="tab10"); ax[1].set_title(f"AFTER ({nrow}x{ncol} grid)"); ax[1].set_aspect("equal"); ax[1].invert_yaxis()
        plt.tight_layout(); plt.savefig(qc,dpi=110); print(f"[qc] {qc}")

if __name__=="__main__":
    ap=argparse.ArgumentParser(); ap.add_argument("path"); ap.add_argument("--out"); ap.add_argument("--ncol",type=int,default=3); ap.add_argument("--qc")
    a=ap.parse_args(); regrid(a.path,a.out,a.ncol,qc=a.qc)
