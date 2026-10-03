#!/usr/bin/env python
# Remove large empty FOV-layout gaps from a sample viewer HTML (cells + mosaic together).
# Bakes the render-time y-flip into the data (global True->False) so x and y are compacted
# in one consistent display space, then crops the matching strips out of the rawImage mosaic.
import json, re, sys, base64, argparse
from io import BytesIO
import numpy as np, pandas as pd
try:
    from PIL import Image
    Image.MAX_IMAGE_PIXELS=None
except Exception:
    Image=None

def extract_data(txt):
    m=re.search(r"const DATA\s*=\s*", txt); s=m.end(); depth=0; i=s
    while i<len(txt):
        c=txt[i]
        if c=='{': depth+=1
        elif c=='}':
            depth-=1
            if depth==0: i+=1; break
        i+=1
    return s,i,json.loads(txt[s:i])

def empty_bands(coords, binsize=50, min_gap=500):
    lo,hi=float(coords.min()),float(coords.max())
    edges=np.arange(lo,hi+binsize,binsize)
    h,_=np.histogram(coords,bins=edges); empty=h==0; bands=[]; i=0
    while i<len(empty):
        if empty[i]:
            j=i
            while j<len(empty) and empty[j]: j+=1
            if (j-i)*binsize>=min_gap: bands.append((float(edges[i]),float(edges[j])))
            i=j
        else: i+=1
    return bands

def typical_pitch(coords, min_gap=500):
    """median FOV-centroid spacing below min_gap (the 'normal' column distance)."""
    s=np.sort(np.array(sorted(set(np.round(coords/50)*50))))
    cols=[]
    for v in s:
        if not cols or v-cols[-1]>150: cols.append(v)
    d=np.diff(cols); small=d[d<min_gap]
    return float(np.median(small)) if len(small) else 400.0

def compact_axis(coords, bands, keep):
    """shift cells to shrink each empty band to `keep`; return new coords + removed world sub-ranges."""
    out=coords.astype(float).copy(); removed=[]
    for a,b in sorted(bands):
        rem=(b-a)-keep
        if rem<=0: continue
        out[coords>b]-=rem
        removed.append((a+keep,b))   # world range physically removed
    return out, removed

def crop_image(datauri, x_removed, y_removed, ext):
    if Image is None or not datauri: return datauri, ext
    b64=datauri.split(",",1)[1]
    img=Image.open(BytesIO(base64.b64decode(b64))).convert("RGB"); W,H=img.size
    e0,e1,e2,e3=ext
    def kept_px(removed, lo, hi, npix):
        rm=sorted([((a-lo)/(hi-lo)*npix,(b-lo)/(hi-lo)*npix) for a,b in removed])
        kept=[]; cur=0
        for a,b in rm:
            a=max(0,int(round(a))); b=min(npix,int(round(b)))
            if a>cur: kept.append((cur,a))
            cur=max(cur,b)
        if cur<npix: kept.append((cur,npix))
        return kept
    xk=kept_px(x_removed,e0,e2,W); yk=kept_px(y_removed,e1,e3,H)
    strips=[img.crop((a,0,b,H)) for a,b in xk]; nW=sum(s.width for s in strips)
    ci=Image.new("RGB",(max(1,nW),H)); x=0
    for s in strips: ci.paste(s,(x,0)); x+=s.width
    strips=[ci.crop((0,a,ci.width,b)) for a,b in yk]; nH=sum(s.height for s in strips)
    fi=Image.new("RGB",(ci.width,max(1,nH))); y=0
    for s in strips: fi.paste(s,(0,y)); y+=s.height
    e2n=e2-sum(b-a for a,b in x_removed); e3n=e3-sum(b-a for a,b in y_removed)
    buf=BytesIO(); fi.save(buf,format="JPEG",quality=84)
    return "data:image/jpeg;base64,"+base64.b64encode(buf.getvalue()).decode(), [e0,e1,e2n,e3n]

def process(path, out=None, min_gap=500, qc=None):
    txt=open(path,encoding="utf-8").read(); s,e,D=extract_data(txt)
    x=np.array(D["x"],float); y=np.array(D["y"],float)
    x0,y0=x.copy(),y.copy(); fov=np.array(D.get("fov",[0]*len(x)))
    # bake render-time y-flip so cells + rawExtent share one display space
    if D.get("global"):
        a,b=y.min(),y.max(); y=b-(y-a); D["global"]=False
    kx,ky=typical_pitch(x,min_gap),typical_pitch(y,min_gap)
    xb,yb=empty_bands(x,min_gap=min_gap),empty_bands(y,min_gap=min_gap)
    xn,xr=compact_axis(x,xb,kx); yn,yr=compact_axis(y,yb,ky)
    if "rawImage" in D and D["rawImage"]:
        D["rawImage"],D["rawExtent"]=crop_image(D["rawImage"],xr,yr,D["rawExtent"])
    D["x"]=[round(float(v),1) for v in xn]; D["y"]=[round(float(v),1) for v in yn]
    newjson=json.dumps(D,separators=(",",":"))
    out=out or path
    open(out,"w",encoding="utf-8").write(txt[:s]+newjson+txt[e:])
    print(f"[degap] {path.split('/')[-1]}: x-gaps {len(xr)} (removed {sum(b-a for a,b in xr):.0f}um) "
          f"y-gaps {len(yr)} (removed {sum(b-a for a,b in yr):.0f}um) keepX {kx:.0f} keepY {ky:.0f}")
    if qc:
        import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
        yy=(y0.max()-(y0-y0.min())) if False else y  # display y already computed
        fig,ax=plt.subplots(1,2,figsize=(13,5))
        ax[0].scatter(x0,(y0.max()+y0.min())-y0 if D.get("global") is False else y0,s=1,c=pd.factorize(fov)[0],cmap="tab20"); ax[0].set_title("BEFORE (empty gaps)"); ax[0].set_aspect("equal")
        ax[1].scatter(xn,yn,s=1,c=pd.factorize(fov)[0],cmap="tab20"); ax[1].set_title("AFTER (gaps compacted)"); ax[1].set_aspect("equal")
        plt.tight_layout(); plt.savefig(qc,dpi=110); print(f"[qc] {qc}")

if __name__=="__main__":
    ap=argparse.ArgumentParser(); ap.add_argument("path"); ap.add_argument("--out"); ap.add_argument("--min-gap",type=float,default=500); ap.add_argument("--qc")
    a=ap.parse_args(); process(a.path,a.out,a.min_gap,a.qc)
