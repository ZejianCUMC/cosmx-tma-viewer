#!/usr/bin/env python
# =============================================================================
# Build sample_viewer_template_cores.html from the golden template
# =============================================================================
# Description : Adds the per-CORE annotation layer to the viewer shell:
#               (1) a frame per TMA core, drawn clear of every dot, with the
#                   core_id inline in the top border (----<core_id>----);
#               (2) frame colour by treatment timepoint (pre = azure,
#                   post = deep red) chosen to not collide with the L1
#                   Epithelial (#e6194B) / Immune (#4363d8) dot colours;
#               (3) a per-core tag row in the header carrying the core's own
#                   clinical string + a show/hide checkbox (default shown);
#               (4) hover a core label -> clinical tooltip.
#               Every edit is anchored on an exact substring of the original
#               template and asserted to apply exactly once, so the shell stays
#               byte-compatible with the build/degap/regrid DATA-span contract.
# Input       : sample_viewer_template.html  (golden shell, must contain __DATA__)
# Output      : sample_viewer_template_cores.html
# Conda env   : any python3
# Validated   : 2026-09-15 (<unit> test viewer)
# Parameters  : See DEFAULTS dict below
# =============================================================================

import argparse, os, re

HERE = os.path.dirname(os.path.abspath(__file__))
DEFAULTS = {
    "src": os.path.join(HERE, "sample_viewer_template.html"),   # golden shell
    "out": os.path.join(HERE, "sample_viewer_template_cores.html"),
    # frame colours: deliberately NOT the L1 dot colours (#e6194B / #4363d8)
    "col_pre":  "#00C2FF",   # azure — vs L1 Immune indigo #4363d8
    "col_post": "#B00000",   # deep red — vs L1 Epithelial crimson #e6194B
}

# ---- (1) CSS for the core tag row -------------------------------------------
CSS_ANCHOR = ".note{color:var(--mut);font-size:10.5px;margin-top:3px}"
CSS_ADD = """
.corerow{display:flex;gap:6px;flex-wrap:wrap;align-items:center}
.crhdr{color:var(--mut);font-size:11px}
.ctag{display:inline-flex;align-items:center;gap:6px;background:var(--chip);border:1px solid var(--edge);border-left:4px solid var(--cc,var(--mut));border-radius:6px;padding:3px 8px;font-size:11px;white-space:nowrap;cursor:pointer;user-select:none}
.ctag:hover{border-color:var(--accent);border-left-color:var(--cc,var(--mut))}
.ctag b{color:var(--cc,var(--fg));font-weight:700}
.ctag .cclin{color:var(--fg)}
.ctag .cnt{color:var(--mut);font-size:10px}
.ctag.off{opacity:.42}
.ctag.off .cclin{text-decoration:line-through}
.ctag.na{opacity:.55;cursor:default;font-style:italic}
.ctag input{margin:0;cursor:pointer}
.crlegend{display:inline-flex;align-items:center;gap:5px;color:var(--mut);font-size:10.5px}
.crlegend i{width:15px;border-top:3px solid;display:inline-block;vertical-align:middle}
.crbtn{font-size:10.5px;padding:2px 7px}"""

# ---- (2) header DOM: frame toggle + core tag row -----------------------------
HDR_ANCHOR = '  <button id="reset">Reset view</button>\n'
HDR_ADD = ('  <span class="gc"><label style="cursor:pointer;user-select:none">'
           '<input type="checkbox" id="framechk" checked style="vertical-align:middle">'
           ' Core frames</label></span>\n')
PATROW_ANCHOR = ' <div class="patrow" id="patinfo"></div>\n'
PATROW_ADD = ' <div class="corerow" id="coreinfo"></div>\n'

# ---- (3) JS: core model, frame drawing, tags, tooltip ------------------------
JS_ANCHOR = "// patient header\n"
JS_ADD = """// --- TMA core frames + per-core clinical annotation ---------------------------
// One physical tissue piece = one frame. Sub-cores of the same punch SHARE the
// core_id, so they share one label, one clinical record, one header tag and one
// show/hide checkbox (<unit> <core_id> is two pieces ~3mm apart on the slide).
const CORES=DATA.cores||[];
const CORE_TPCOL={pre_treatment:"__COL_PRE__",post_treatment:"__COL_POST__"};
const coreCol=c=>CORE_TPCOL[c.tp]||"#9aa4b2";
const CORE_IDS=[...new Set(CORES.map(c=>c.id))];
const shownCore={};CORE_IDS.forEach(id=>shownCore[id]=true);
let showFrames=true;
const COREIDX=new Int16Array(N);
(function(){const mp=DATA.coreFov||{};for(let i=0;i<N;i++){const v=mp[DATA.fov[i]];COREIDX[i]=(v===undefined?-1:v)}})();
const coreHidden=i=>{const c=COREIDX[i];return c>=0&&!shownCore[CORES[c].id]};
function coreAgg(id){const ss=CORES.filter(c=>c.id===id);return{c:ss[0],n:ss.reduce((a,b)=>a+b.n,0),nfov:ss.reduce((a,b)=>a+b.nfov,0),parts:ss.length}}
function coreTip(c){const g=coreAgg(c.id);const pc=c.nparts>1?` <span style="color:var(--mut)">piece ${c.part} of ${c.nparts}</span>`:"";
 return `<b style="color:${coreCol(c)}">${c.id}</b>${pc}<br>Path: ${c.path_dx}<br>${c.tpLabel}<br>Stage ${c.stage} \u00b7 Grade ${c.grade}<br><span style="color:var(--mut)">Block ${c.block} \u00b7 ${c.n.toLocaleString()} cells \u00b7 ${c.nfov} FOVs${c.nparts>1?` \u00b7 core total ${g.n.toLocaleString()}`:""}</span>`}
function maskHiddenCores(p,bg){if(!CORES.length)return;const ctx=p.ctx;ctx.save();ctx.setTransform(DPR,0,0,DPR,0,0);ctx.fillStyle=bg;CORES.forEach(c=>{if(shownCore[c.id])return;const s=c.slot;ctx.fillRect(SX(s[0]),SY(s[1]),(s[2]-s[0])*view.s,(s[3]-s[1])*view.s)});ctx.restore()}
function drawCores(p){p.labelRects=[];if(!CORES.length||!showFrames)return;const ctx=p.ctx;ctx.save();ctx.setTransform(DPR,0,0,DPR,0,0);
 const grow=PSZ+2;ctx.font='bold 12px -apple-system,BlinkMacSystemFont,"Segoe UI",Roboto,sans-serif';ctx.textBaseline="middle";ctx.textAlign="center";
 CORES.forEach((c,ci)=>{if(!shownCore[c.id])return;const col=coreCol(c);
  const X0=SX(c.box[0])-grow,Y0=SY(c.box[1])-grow,X1=SX(c.box[2])+grow,Y1=SY(c.box[3])+grow;
  if(X1<-60||X0>p.W+60||Y1<-60||Y0>p.H+60)return;
  const tw=ctx.measureText(c.id).width,half=tw/2+7;
  let lx=(X0+X1)/2;const loA=Math.max(X0,0)+half+2,hiA=Math.min(X1,p.W)-half-2;
  if(hiA>loA)lx=Math.min(Math.max(lx,loA),hiA);
  const ly=Math.min(Math.max(Y0,15),p.H-7),clamped=Math.abs(ly-Y0)>0.5;
  ctx.strokeStyle=col;ctx.lineWidth=2;ctx.setLineDash([]);
  ctx.beginPath();                                  // 3 full sides + broken top
  ctx.moveTo(X0,Y0);ctx.lineTo(X0,Y1);ctx.lineTo(X1,Y1);ctx.lineTo(X1,Y0);
  if(clamped){ctx.lineTo(X0,Y0)}else{ctx.lineTo(lx+half,Y0);ctx.moveTo(lx-half,Y0);ctx.lineTo(X0,Y0)}
  ctx.stroke();
  const bw=tw+12,bh=15;                             // label pill (legible over morphology)
  ctx.fillStyle=getComputedStyle(document.body).backgroundColor||"#0f1216";
  ctx.fillRect(lx-bw/2,ly-bh/2,bw,bh);
  if(clamped){ctx.strokeStyle=col;ctx.lineWidth=1;ctx.strokeRect(lx-bw/2,ly-bh/2,bw,bh)}
  ctx.fillStyle=col;ctx.fillText(c.id,lx,ly+0.5);
  p.labelRects.push({x0:lx-bw/2,y0:ly-bh/2,x1:lx+bw/2,y1:ly+bh/2,ci:ci})});
 ctx.restore()}
function coreLabelAt(p,mx,my){const L=p.labelRects||[];for(let i=L.length-1;i>=0;i--){const r=L[i];if(mx>=r.x0&&mx<=r.x1&&my>=r.y0&&my<=r.y1)return r.ci}return -1}
function renderCoreTags(){const box=document.getElementById("coreinfo");if(!box)return;
 if(!CORES.length){box.style.display="none";return}
 const legend=`<span class="crlegend"><i style="border-color:${CORE_TPCOL.pre_treatment}"></i>Pre-treatment<i style="border-color:${CORE_TPCOL.post_treatment};margin-left:6px"></i>Post-treatment</span>`;
 const tags=CORE_IDS.map((id,k)=>{const g=coreAgg(id),c=g.c,on=shownCore[id];
  const pc=g.parts>1?` \u00b7 ${g.parts} pieces`:"";
  return `<label class="ctag${on?"":" off"}" style="--cc:${coreCol(c)}" title="${id} \u2014 ${c.clin} \u2014 block ${c.block}${pc}"><input type="checkbox" data-ck="${k}"${on?" checked":""}><b>${id}</b><span class="cclin">(${c.clin})</span><span class="cnt">${g.n.toLocaleString()} cells \u00b7 ${g.nfov} FOVs${pc}</span></label>`}).join("");
 box.innerHTML=`<span class="crhdr">Cores (${CORE_IDS.length})</span>${tags}<button class="crbtn" id="corall">All</button><button class="crbtn" id="cornone">None</button>${legend}`;
 box.querySelectorAll("input[data-ck]").forEach(inp=>inp.onchange=()=>{shownCore[CORE_IDS[+inp.dataset.ck]]=inp.checked;renderCoreTags();drawAll()});
 const all=document.getElementById("corall"),non=document.getElementById("cornone");
 if(all)all.onclick=()=>{CORE_IDS.forEach(id=>shownCore[id]=true);renderCoreTags();drawAll()};
 if(non)non.onclick=()=>{CORE_IDS.forEach(id=>shownCore[id]=false);renderCoreTags();drawAll()}}
"""

# ---- (4) hide hidden-core cells in the scatter loop ---------------------------
LOOP_ANCHOR = ('const _hi=[];for(let i=0;i<N;i++){if(f.kind==="cat"&&p.state)'
               '{const _st=p.state.get(vals[i])||0;')
LOOP_NEW = ('const _hi=[];for(let i=0;i<N;i++){if(coreHidden(i))continue;'
            'if(f.kind==="cat"&&p.state){const _st=p.state.get(vals[i])||0;')

# ---- (5) draw frames in every render path ------------------------------------
RAWIMG_ANCHOR = ('(e[2]-e[0])*view.s,(e[3]-e[1])*view.s)}drawRegion(p);legend(p,f);return}')
RAWIMG_NEW = ('(e[2]-e[0])*view.s,(e[3]-e[1])*view.s)}maskHiddenCores(p,"#000");'
              'drawRegion(p);drawCores(p);legend(p,f);return}')
END_ANCHOR = 'ctx.globalAlpha=1;drawRegion(p);legend(p,f,rng)}'
END_NEW = 'ctx.globalAlpha=1;drawRegion(p);drawCores(p);legend(p,f,rng)}'

# ---- (6) header info line: cores instead of one core's path_dx ----------------
FINFO_ANCHOR = ('document.getElementById("finfo").textContent=`${DATA.sample||DATA.fov} · '
                '${N.toLocaleString()} cells${DATA.fovs?" · "+DATA.fovs+" FOVs":""} · '
                '${(DATA.meta.pathdx&&DATA.meta.pathdx[0])||""}`;')
FINFO_NEW = ('document.getElementById("finfo").textContent=`${DATA.sample||DATA.fov} · '
             '${N.toLocaleString()} cells${DATA.fovs?" · "+DATA.fovs+" FOVs":""}'
             '${CORE_IDS.length?" · "+CORE_IDS.length+" cores"+(CORES.length>CORE_IDS.length?'
             '" ("+CORES.length+" pieces)":""):""}`;renderCoreTags();')

# ---- (7) hover: core label tooltip first, then cells (skipping hidden) --------
HOVER_ANCHOR = ('function hover(e){const p=panels.find(q=>q.cv===e.target);if(!p)return;'
                'const r=p.cv.getBoundingClientRect(),mx=e.clientX-r.left,my=e.clientY-r.top;'
                'let best=-1,bd=64;for(let i=0;i<N;i++){')
HOVER_NEW = ('function hover(e){const p=panels.find(q=>q.cv===e.target);if(!p)return;'
             'const r=p.cv.getBoundingClientRect(),mx=e.clientX-r.left,my=e.clientY-r.top;'
             'const _lc=coreLabelAt(p,mx,my);p.cv.style.cursor=_lc>=0?"pointer":"crosshair";'
             'if(_lc>=0){tt.innerHTML=coreTip(CORES[_lc]);tt.style.left=(e.clientX+14)+"px";'
             'tt.style.top=(e.clientY+14)+"px";tt.style.display="block";return}'
             'let best=-1,bd=64;for(let i=0;i<N;i++){if(coreHidden(i))continue;')

# ---- (8) cell tooltip: name the core the cell belongs to ---------------------
TIP_ANCHOR = 'tt.innerHTML=`<b>${DATA.type[best]}</b> / ${DATA.subtype[best]}<br>'
TIP_NEW = ('tt.innerHTML=`${COREIDX[best]>=0?`<span style="color:${coreCol(CORES[COREIDX[best]])}">'
           '${CORES[COREIDX[best]].id}</span> · ${CORES[COREIDX[best]].path_dx} · '
           '${CORES[COREIDX[best]].tpLabel}<br>`:""}<b>${DATA.type[best]}</b> / ${DATA.subtype[best]}<br>')

# ---- (9) frame toggle handler ------------------------------------------------
TAIL_ANCHOR = 'window.addEventListener("resize",resize);resize();'
TAIL_NEW = ('(function(){const c=document.getElementById("framechk");if(!c)return;'
            'if(!CORES.length){c.disabled=true;c.parentElement.style.opacity=0.45;'
            'c.parentElement.title="No core annotation in this viewer";showFrames=false;return}'
            'c.onchange=e=>{showFrames=e.target.checked;drawAll()}})();\n'
            'window.addEventListener("resize",resize);resize();')


# ---- (10) patient chip row: aggregate Grade / Stage / Treatment over the CORES
#      actually present (post-QC), each with its core count -----------------------
PAT_LINE_PREFIX = 'document.getElementById("patinfo").innerHTML='
PAT_NEW = ('const _nc=n=>`${n} core${n>1?"s":""}`;\n'
           'function coreSummary(){const u=CORE_IDS.map(id=>CORES.find(c=>c.id===id));\n'
           ' const agg=(fn,order)=>{const m={};u.forEach(c=>{const v=fn(c);if(v)m[v]=(m[v]||0)+1});\n'
           '  const ks=Object.keys(m).sort((a,b)=>(order?order.indexOf(a)-order.indexOf(b):0)||m[b]-m[a]||a.localeCompare(b));\n'
           '  return ks.map(k=>`${k}(${_nc(m[k])})`).join(", ")};\n'
           ' return [["Grade",agg(c=>c.grade)],["Stage",agg(c=>c.stage)],\n'
           '  ["Treatment",agg(c=>c.tp==="pre_treatment"?"Pre":c.tp==="post_treatment"?"Post":"",'
           '["Pre","Post"])]]}\n'
           'function renderPatChips(){const chip=(k,v,w)=>`<span class="chip${w?" warn":""}">'
           '<b>${k}</b>${v}</span>`;const out=[];let done=false;\n'
           ' const sum=()=>{done=true;coreSummary().forEach(([k,v])=>{if(v)out.push(chip(k,v))})};\n'
           ' Object.entries(PAT).filter(([k,v])=>v&&v!=="NA"&&v!=="nan").forEach(([k,v])=>{\n'
           '  if(k==="Age")v=parseFloat(v).toFixed(0)+"y";\n'
           '  if(k==="TTR (months)")v=parseFloat(v).toFixed(1)+" mo";\n'
           '  const w=(k==="Progression"&&v==="Yes")||(k==="Outcome"&&/PROGRE|DECEAS/i.test(v));\n'
           '  out.push(chip(k,v,w));if(k==="Age")sum()});\n'
           ' if(!done)sum();document.getElementById("patinfo").innerHTML=out.join("")}\n'
           'renderPatChips();')

# ---- (11) hide the important region control entirely when this unit has no important region --------
REGION_ANCHOR = ('(function(){const _c=document.getElementById("regionchk");if(!_c)return;'
              'if(!DATA.importantRegion){_c.disabled=true;_c.parentElement.style.opacity=0.45;'
              '_c.parentElement.title="No important region for this patient"}'
              '_c.onchange=e=>{showRegion=e.target.checked;drawAll()}})();')
REGION_NEW = ('(function(){const _c=document.getElementById("regionchk");if(!_c)return;'
           'if(!DATA.importantRegion){const _w=_c.closest(".gc")||_c.parentElement;'
           '_w.style.display="none";return}'
           '_c.onchange=e=>{showRegion=e.target.checked;drawAll()}})();')

EDITS = [("CSS", CSS_ANCHOR, CSS_ANCHOR + CSS_ADD),
         ("header frame toggle", HDR_ANCHOR, HDR_ADD + HDR_ANCHOR),
         ("core tag row", PATROW_ANCHOR, PATROW_ANCHOR + PATROW_ADD),
         ("core JS block", JS_ANCHOR, JS_ADD + JS_ANCHOR),
         ("scatter skip hidden", LOOP_ANCHOR, LOOP_NEW),
         ("rawimg frames", RAWIMG_ANCHOR, RAWIMG_NEW),
         ("scatter frames", END_ANCHOR, END_NEW),
         ("finfo line", FINFO_ANCHOR, FINFO_NEW),
         ("hover core label", HOVER_ANCHOR, HOVER_NEW),
         ("cell tooltip core", TIP_ANCHOR, TIP_NEW),
         ("frame toggle handler", TAIL_ANCHOR, TAIL_NEW),
         ("hide empty important region control", REGION_ANCHOR, REGION_NEW)]


def main():
    ap = argparse.ArgumentParser(description="Build the core-annotated viewer template.")
    for k, v in DEFAULTS.items():
        ap.add_argument("--" + k.replace("_", "-"), default=v)
    a = ap.parse_args()

    t = open(a.src, encoding="utf-8").read()
    assert "__DATA__" in t, "source template lost its __DATA__ placeholder"

    # whole-line edit: the patinfo builder is retyped nowhere, it is located in
    # the source and replaced end-to-end (its exact minified text must not be
    # transcribed by hand -- a single space difference silently breaks matching)
    lines = t.split("\n")
    hits = [i for i, ln in enumerate(lines) if ln.startswith(PAT_LINE_PREFIX)]
    assert len(hits) == 1, f"patinfo line matched {len(hits)} times (expected 1)"
    lines[hits[0]] = PAT_NEW
    t = "\n".join(lines)

    for name, old, new in EDITS:
        n = t.count(old)
        assert n == 1, f"anchor '{name}' matched {n} times (expected 1)"
        t = t.replace(old, new, 1)
    t = t.replace("__COL_PRE__", a.col_pre).replace("__COL_POST__", a.col_post)

    assert t.count("__DATA__") == 1
    assert "__COL_" not in t
    # the DATA-span contract every postprocess script depends on
    assert re.search(r"const DATA\s*=\s*__DATA__", t), "DATA assignment shape changed"
    # feature-checker literals that must survive (verify_viewer_features.py)
    for lit in ['id="psz"', 'id="pop"', 'id="reset"', 'id="ci0"', 'id="ci1"', 'id="cp0"',
                'id="cp1"', '.li[data-cat]', 'p.state', 'DATA.rawImage?"__rawimg__"',
                'function drawRegion(p)', 'showRegion', 'id="regionchk"', 'Important region',
                'CellComposite morphology']:
        assert lit in t, f"feature-checker literal lost: {lit}"
    open(a.out, "w", encoding="utf-8").write(t)
    print(f"[template] {a.out}  ({len(t)} bytes, +{len(t)-len(open(a.src,encoding='utf-8').read())} "
          f"vs golden)  pre={a.col_pre} post={a.col_post}")
    print(f"[template] {len(EDITS)} anchored edits applied, all feature-checker literals intact")


if __name__ == "__main__":
    main()
