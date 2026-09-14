"""
Match each simplified solid to a full-assembly solid by bounding-box overlap.

    "C:/Program Files/FreeCAD 1.1/bin/python.exe" scripts/investigation/06_match_simplified_to_full.py

Writes output/simp_to_full.json: {simp, full, score} per simplified solid, where
score is intersection volume over the larger of the two boxes (1.0 = identical
boxes, 0 = no overlap).

This is the hop that carries material across files: the simplified STEP has no
material data at all, so each of its solids has to be identified with a solid in
the colour-coded full assembly. Both files are in the same CAD frame, so raw
overlap is meaningful without any transform.

Deliberately kernel-only -- no parsed coordinates anywhere -- so it is
independent of the text-parsing path in 01. If the two disagreed, that would be
a real signal rather than a shared bug.

The quality histogram it prints (>=0.99 / >=0.90 / >=0.50 / <0.50) is the
honest measure of the match, and the poor matches are listed. Two solids match
badly because simplification changed them substantially; the extractor carries a
material_match score per solid for exactly this reason.

Slow: reads both files with the kernel, the full one being 74MB.
"""

import collections
import json
import os
import time

import FreeCAD  # noqa: F401
import Part

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
OUTDIR = os.path.join(ROOT, "output")
SIMPLE = os.path.join(ROOT, "STM_STP_files",
                      "F10269585--_1-G4 Shield House Simplified.stp")
FULL = os.path.join(ROOT, "STM_STP_files", "F10258491--_1-Shield House Square.stp")

t = time.time()
S = Part.read(SIMPLE)
F = Part.read(FULL)
print('read sec', round(time.time() - t, 1), 'simp', len(S.Solids), 'full', len(F.Solids))

fb = []
for j, s in enumerate(F.Solids):
    b = s.BoundBox
    fb.append((j, b.XMin, b.XMax, b.YMin, b.YMax, b.ZMin, b.ZMax, s.Volume))


def ov(a0, a1, b0, b1):
    return max(0.0, min(a1, b1) - max(a0, b0))


res = []
for i, s in enumerate(S.Solids):
    b = s.BoundBox
    sv = b.XLength * b.YLength * b.ZLength
    best = (0.0, None)
    for (j, x0, x1, y0, y1, z0, z1, fvol) in fb:
        iv = (ov(b.XMin, b.XMax, x0, x1) * ov(b.YMin, b.YMax, y0, y1)
              * ov(b.ZMin, b.ZMax, z0, z1))
        if iv <= 0:
            continue
        fv = (x1 - x0) * (y1 - y0) * (z1 - z0)
        sc = iv / max(sv, fv, 1e-9)
        if sc > best[0]:
            best = (sc, j)
    res.append((i, best[1], best[0]))

q = collections.Counter()
for i, j, sc in res:
    q['>=0.99' if sc >= 0.99 else '>=0.90' if sc >= 0.90
      else '>=0.50' if sc >= 0.5 else '<0.50'] += 1
print('=== kernel-only match quality ===')
for k in ('>=0.99', '>=0.90', '>=0.50', '<0.50'):
    print('   %-8s %4d' % (k, q[k]))

json.dump([{'simp': i, 'full': j, 'score': round(sc, 4)} for i, j, sc in res],
          open(os.path.join(OUTDIR, 'simp_to_full.json'), 'w'), indent=1)
print('wrote output/simp_to_full.json')
print('poor matches:', [(i, round(sc, 3)) for i, j, sc in res if sc < 0.5][:15])
