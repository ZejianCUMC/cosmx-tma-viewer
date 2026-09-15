#!/usr/bin/env python
"""Verify a re-stitched viewer changed ONLY the morphology vs a reference viewer.
Usage: python verify_only_morphology_changed.py NEW.html REF.html"""
import sys,re,json,base64,io,hashlib
from PIL import Image
Image.MAX_IMAGE_PIXELS=None
def data(p):
    s=open(p,encoding="utf-8",errors="ignore").read();i=s.find('DATA=');j=s.find('{',i);d=0;k=j
    while k<len(s):
        c=s[k]
        if c=='{':d+=1
        elif c=='}':
            d-=1
            if d==0:break
        k+=1
    return json.loads(s[j:k+1])
def rd(D):
    raw=base64.b64decode(D["rawImage"].split(",",1)[1]);im=Image.open(io.BytesIO(raw));return im.size,len(raw)
def h(x):return hashlib.md5(json.dumps(x,sort_keys=True,separators=(',',':')).encode()).hexdigest()[:8]
N,R=data(sys.argv[1]),data(sys.argv[2])
(nw,nh),nb=rd(N);(rw,rh),rb=rd(R)
print(f"morphology: new {nw}x{nh} ({nb//1024//1024}MB) vs ref {rw}x{rh} ({rb//1024}KB)")
print("rawExtent equal:",[round(v,1) for v in N['rawExtent']]==[round(v,1) for v in R['rawExtent']])
ok=all(h(N.get(k))==h(R.get(k)) for k in ["subL2","subtype","type","x","y","fov","tumor","genes","channels","meta"])
print("ALL non-morphology fields identical:",ok)
sys.exit(0 if ok else 1)
