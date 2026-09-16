#!/usr/bin/env python
"""Splice a built viewer's DATA span into the core-annotated shell.
Usage: splice.py <built.html> <cores_template.html> <out.html>"""
import io, re, sys

def data_span(txt):
    m = re.search(r"const DATA\s*=\s*", txt)
    s = m.end(); d = 0; i = s
    while i < len(txt):
        c = txt[i]
        if c == '{': d += 1
        elif c == '}':
            d -= 1
            if d == 0: i += 1; break
        i += 1
    return s, i

src, tpl_p, out = sys.argv[1:4]
txt = io.open(src, encoding="utf-8").read()
s, e = data_span(txt)
tpl = io.open(tpl_p, encoding="utf-8").read()
assert tpl.count("__DATA__") == 1, "template lost its __DATA__ placeholder"
io.open(out, "w", encoding="utf-8").write(tpl.replace("__DATA__", txt[s:e]))
print(f"[splice] {out.split('/')[-1]}  {(len(tpl)+e-s)/1e6:.2f} MB")
