#!/usr/bin/env python
# =============================================================================
# cosmx-tma-viewer : headless-browser verification of built viewers
# =============================================================================
# Description : Opens each viewer in headless Chrome (DevTools protocol) and checks what a user actually sees:
#               no JS errors; canvas drawing size == displayed size; real mouse at the ON-SCREEN position of cells
#               finds the right cell (zoomed in, largest dot size); hidden-by-default categories are neither drawn
#               nor hovered; no cell tooltip on morphology / H&E layers; grouped level-3 legend complete and inside
#               its panel; region labels inside the canvas and clear of core labels; full-resolution morphology
#               loads only when zoomed in; every cell lies inside its core frame (data check).
# Input       : viewer HTML files; a Chrome / Chromium binary
# Output      : printed table + JSON report; exit code 1 if any check fails
# Key deps    : requests, websocket-client (pip install requests websocket-client)
# =============================================================================
import argparse, json, os, re, shutil, subprocess, sys, tempfile, time
from pathlib import Path
import requests, websocket

DEFAULTS = {"port": 9238, "width": 1700, "height": 1050, "hover_cells": 40, "hover_min_frac": 0.95, "load_timeout_s": 180}

class Page:
    def __init__(self, port):
        pg = [p for p in requests.get(f"http://localhost:{port}/json").json() if p.get("type") == "page"][0]
        self.ws = websocket.create_connection(pg["webSocketDebuggerUrl"], timeout=300); self.seq = 0
    def call(self, m, p=None):
        self.seq += 1; self.ws.send(json.dumps({"id": self.seq, "method": m, "params": p or {}}))
        while True:
            r = json.loads(self.ws.recv())
            if r.get("id") == self.seq: return r.get("result", {})
    def ev(self, e): return self.call("Runtime.evaluate", {"expression": e, "returnByValue": True})["result"].get("value")
    def mouse(self, x, y): self.call("Input.dispatchMouseEvent", {"type": "mouseMoved", "x": x, "y": y})

def data_check(fn):
    t = Path(fn).read_text(encoding="utf-8"); m = re.search(r"const DATA\s*=\s*", t); D = json.JSONDecoder().raw_decode(t, m.end())[0]
    bad = 0
    for i, (x, y, f) in enumerate(zip(D["x"], D["y"], D["fov"])):
        c = D["cores"][D["coreFov"][str(f)]]; b = c["box"]
        bad += not (b[0] - 1e-6 <= x <= b[2] + 1e-6 and b[1] - 1e-6 <= y <= b[3] + 1e-6)
    return bad

def check(page, fn, P):
    R = {"viewer": Path(fn).name}
    page.call("Page.navigate", {"url": Path(fn).resolve().as_uri()})
    t0 = time.time()
    while time.time() - t0 < P["load_timeout_s"]:
        if page.ev("typeof panels!=='undefined'&&document.readyState==='complete'&&(!RAWIMG||RAWIMG.complete)"): break
        time.sleep(0.2)
    time.sleep(1.5)
    R["js_errors"] = page.ev("window.__errs||[]")
    R["canvas_sync"] = page.ev("panels.every(p=>{const r=p.cv.getBoundingClientRect();return Math.abs(r.width-p.W)<1&&Math.abs(r.height-p.H)<1})")
    # level-3 view on the right panel, hidden-by-default state, count drawn dots
    page.ev("(()=>{const p=panels[1];p.feat='__subtype__';p.state=hiddenState();drawAll()})()")
    page.ev("window._arc=0;const _o=panels[1].ctx.arc.bind(panels[1].ctx);panels[1].ctx.arc=function(){window._arc++;return _o(...arguments)}")
    exp = page.ev("(()=>{let n=0;for(let i=0;i<N;i++){if(coreHidden(i)||LQC[i])continue;const X=SX(XS[i]),Y=SY(YS[i]);if(X<-5||X>panels[1].W+5||Y<-5||Y>panels[1].H+5)continue;n++}return n})()")
    page.ev("window._arc=0;drawP(panels[1])"); R["hidden_not_drawn"] = page.ev("window._arc") == exp
    # legend: every level-3 value listed once, nothing overflows the legend box
    R["legend_ok"] = page.ev("""(()=>{const lg=document.getElementById('lg1'),L=lg.getBoundingClientRect();const it=[...lg.querySelectorAll('.l2h,.li[data-cat]')];
      const miss=L3U.filter(v=>!lg.querySelector('[data-cat="'+CSS.escape(v)+'"]')&&!lg.querySelector('.l2h[data-l2="'+CSS.escape(v)+'"]')).length;
      return miss===0&&it.every(e=>e.getBoundingClientRect().right<=L.right+0.5)&&lg.scrollWidth<=lg.clientWidth+1})()""")
    # zoom on the first core, largest dots, real mouse at on-screen dot centres
    page.ev("(()=>{const s=document.getElementById('psz');s.value=s.max;s.dispatchEvent(new Event('input'));const c=CORES[0],b=c.box,cx=(b[0]+b[2])/2,cy=(b[1]+b[3])/2;view.s=Math.min(panels[1].W,panels[1].H)/((b[3]-b[1])/2.5);view.tx=panels[1].W/2-cx*view.s;view.ty=panels[1].H/2-cy*view.s;drawAll()})()")
    time.sleep(0.4)
    pts = page.ev(f"""(()=>{{const p=panels[1],b=p.cv.getBoundingClientRect(),kx=b.width/p.W,ky=b.height/p.H,out=[];for(let i=0;i<N;i++){{if(LQC[i]||coreHidden(i))continue;const X=SX(XS[i]),Y=SY(YS[i]);if(X>30&&X<p.W-30&&Y>30&&Y<p.H-30)out.push([b.left+X*kx,b.top+Y*ky,DATA.subtype[i]])}}return out.filter((_,k)=>k%Math.max(1,Math.floor(out.length/{P['hover_cells']}))===0).slice(0,{P['hover_cells']})}})()""") or []
    hit = 0
    for x, y, lab in pts:
        page.mouse(x, y); time.sleep(0.01); t = page.ev("tt.style.display!=='none'?tt.innerText:''") or ""
        hit += lab in t
    R["hover_hit"] = f"{hit}/{len(pts)}"; R["hover_ok"] = (not pts) or hit / len(pts) >= P["hover_min_frac"]
    # morphology / H&E layers: no cell tooltip
    page.ev("(()=>{panels[1].feat=FEATS.some(f=>f.k==='__he__')?'__he__':'__rawimg__';panels[0].feat='__rawimg__';drawAll()})()"); time.sleep(0.3)
    cell_tips = 0
    for x, y, _ in pts[:20]:
        page.mouse(x, y); time.sleep(0.01); t = page.ev("tt.style.display!=='none'?tt.innerText:''") or ""
        cell_tips += any(v in t for v in ("FOV ",))
    R["no_cell_tooltip_on_images"] = cell_tips == 0
    # regions: labels inside the canvas and not over core labels
    R["regions"] = page.ev("TLSR.length")
    R["region_labels_ok"] = page.ev("(()=>{drawAll();return panels.every(p=>(p.tlsRects||[]).every(r=>r[0]>=0&&r[1]>=0&&r[2]<=p.W&&r[3]<=p.H&&!(p.labelRects||[]).some(c=>r[0]<c.x1&&r[2]>c.x0&&r[1]<c.y1&&r[3]>c.y0)))})()")
    # full-resolution layer: nothing loaded at fit view, loaded when zoomed
    R["hires"] = page.ev("RAWHI.length")
    if R["hires"]:
        page.ev("(()=>{RAWHI.forEach(h=>{h.img=null});fit();drawAll()})()"); time.sleep(0.3); at_fit = page.ev("RAWHI.filter(h=>h.img).length")
        page.ev("(()=>{const h=RAWHI[0],e=h.ext,cx=(e[0]+e[2])/2,cy=(e[1]+e[3])/2;view.s=2.5;view.tx=panels[0].W/2-cx*view.s;view.ty=panels[0].H/2-cy*view.s;drawAll()})()")
        for _ in range(100):
            if page.ev("RAWHI.some(h=>h.img&&h.img.naturalWidth>0)"): break
            time.sleep(0.2)
        R["hires_ok"] = at_fit == 0 and page.ev("RAWHI.some(h=>h.img&&h.img.naturalWidth>0)")
    R["cells_outside_frames"] = data_check(fn)
    keys = ["canvas_sync", "hidden_not_drawn", "legend_ok", "hover_ok", "no_cell_tooltip_on_images", "region_labels_ok"] + (["hires_ok"] if R["hires"] else [])
    R["pass"] = not R["js_errors"] and all(R[k] for k in keys) and R["cells_outside_frames"] == 0
    return R

def run_pipeline(viewers, chrome=None, port=None, report=None, **kw):
    P = {**DEFAULTS, **kw, **({"port": port} if port else {})}
    chrome = chrome or shutil.which("google-chrome") or shutil.which("chromium") or shutil.which("chromium-browser")
    if not chrome: sys.exit("no Chrome/Chromium found; pass --chrome")
    prof = tempfile.mkdtemp(prefix="tmaviewer_chrome_")
    proc = subprocess.Popen([chrome, "--headless=new", "--disable-gpu", "--no-sandbox", f"--remote-debugging-port={P['port']}",
                             "--remote-allow-origins=*", f"--user-data-dir={prof}", "about:blank"], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    try:
        for _ in range(50):
            try: requests.get(f"http://localhost:{P['port']}/json", timeout=1); break
            except Exception: time.sleep(0.2)
        page = Page(P["port"]); page.call("Page.enable"); page.call("Network.enable"); page.call("Network.setCacheDisabled", {"cacheDisabled": True})
        page.call("Page.addScriptToEvaluateOnNewDocument", {"source": "window.__errs=[];window.addEventListener('error',e=>window.__errs.push(String(e.message)))"})
        page.call("Emulation.setDeviceMetricsOverride", {"width": P["width"], "height": P["height"], "deviceScaleFactor": 1, "mobile": False})
        res = [check(page, f, P) for f in viewers]
    finally:
        proc.terminate(); shutil.rmtree(prof, ignore_errors=True)
    for r in res: print(("PASS " if r["pass"] else "FAIL ") + json.dumps({k: v for k, v in r.items() if k != "pass"}, ensure_ascii=False))
    if report: json.dump(res, open(report, "w"), indent=1)
    print(f"{sum(r['pass'] for r in res)}/{len(res)} viewers pass"); return {"summary": res}

if __name__ == "__main__":
    ap = argparse.ArgumentParser(description="headless-browser checks for built viewers")
    ap.add_argument("viewers", nargs="+"); ap.add_argument("--chrome"); ap.add_argument("--port", type=int); ap.add_argument("--report")
    a = ap.parse_args(); r = run_pipeline(a.viewers, a.chrome, a.port, a.report)
    sys.exit(0 if all(x["pass"] for x in r["summary"]) else 1)
