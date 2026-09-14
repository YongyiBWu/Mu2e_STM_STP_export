"""
Bind full-file kernel solids to the parsed MSB colours by centroid proximity.

    "C:/Program Files/FreeCAD 1.1/bin/python.exe" scripts/investigation/02_bind_colours_by_centroid.py

Reads  output/full_msb_colours.json   (from 01)
Writes output/full_solids_colour.json: kernel solid index -> colour, with an
accurate kernel bounding box and volume.

The point of this step is to move from parsed-point boxes to kernel boxes while
carrying the colour across. Matching is by centroid: MSB entries are bucketed on
a 0.5mm grid and each kernel solid looks in its own bucket and the 26 around it.

This approach turned out NOT to be good enough -- it leaves ambiguity where
solids share a centroid, which is common for concentric parts. Script 05
replaces it with an ordering-based bind that resolves all 224. Kept because the
exact/near/miss counts it prints are what showed the method was insufficient.

Slow: reads the 74MB full assembly with the kernel (~minutes).
"""

import collections
import json
import os
import time

import FreeCAD  # noqa: F401  (must precede Part; see extract_stm_geometry.py)
import Part

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
OUTDIR = os.path.join(ROOT, "output")
FULL = os.path.join(ROOT, "STM_STP_files", "F10258491--_1-Shield House Square.stp")

full = json.load(open(os.path.join(OUTDIR, 'full_msb_colours.json')))
t = time.time()
sh = Part.read(FULL)
print('read sec', round(time.time() - t, 1), 'solids', len(sh.Solids))


def key(x, y, z, q=0.5):
    return (round(x / q), round(y / q), round(z / q))


idx = collections.defaultdict(list)
for e in full:
    idx[key(e['cx'], e['cy'], e['cz'])].append(e)

rows = []
exact = near = miss = 0
for i, s in enumerate(sh.Solids):
    b = s.BoundBox
    c = b.Center
    k0 = key(c.x, c.y, c.z)
    cand = []
    for dx in (-1, 0, 1):
        for dy in (-1, 0, 1):
            for dz in (-1, 0, 1):
                cand.extend(idx.get((k0[0] + dx, k0[1] + dy, k0[2] + dz), []))
    best = None
    for e in cand:
        d = ((e['cx'] - c.x) ** 2 + (e['cy'] - c.y) ** 2 + (e['cz'] - c.z) ** 2) ** 0.5
        if best is None or d < best[0]:
            best = (d, e)
    if best and best[0] < 0.05:
        exact += 1
        col = best[1]['colour']
    elif best and best[0] < 2.0:
        near += 1
        col = best[1]['colour']
    else:
        miss += 1
        col = None
    rows.append({'i': i, 'colour': col,
                 'cx': round(c.x, 3), 'cy': round(c.y, 3), 'cz': round(c.z, 3),
                 'dx': round(b.XLength, 3), 'dy': round(b.YLength, 3),
                 'dz': round(b.ZLength, 3), 'vol': round(s.Volume, 1)})

print('centroid bind: exact(<0.05mm) %d  near(<2mm) %d  miss %d' % (exact, near, miss))
json.dump(rows, open(os.path.join(OUTDIR, 'full_solids_colour.json'), 'w'), indent=1)
print('wrote output/full_solids_colour.json')
print('colour counts:', dict(collections.Counter(r['colour'] for r in rows)))
