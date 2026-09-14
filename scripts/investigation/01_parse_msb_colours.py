"""
Parse the STEP styling chain of the full assembly to get a colour per solid.

    "C:/Program Files/FreeCAD 1.1/bin/python.exe" scripts/investigation/01_parse_msb_colours.py

Writes output/full_msb_colours.json: one entry per MANIFOLD_SOLID_BREP with its
colour name and a bounding box.

Why parse the text rather than ask the kernel: FreeCAD's Part.read gives solids
but drops the STEP presentation layer, which is where colour lives. The chain is

    STYLED_ITEM -> PRESENTATION_STYLE_ASSIGNMENT -> ... -> COLOUR_RGB

with several link types in the middle, so it is walked breadth-first rather than
matched by a fixed pattern. The bounding box here comes from parsed
CARTESIAN_POINTs, not the kernel, so it is only good enough for coarse matching;
02 replaces it with accurate kernel boxes.

This is the first hop of the material recovery described in
extract_stm_geometry.py. Pure text parsing: does not import FreeCAD.
"""

import collections
import json
import os
import re

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
OUTDIR = os.path.join(ROOT, "output")
FULL = os.path.join(ROOT, "STM_STP_files", "F10258491--_1-Shield House Square.stp")

txt = open(FULL, "r", errors="replace").read()
ent = dict(re.findall(r'#(\d+)=([A-Z_0-9]+\([^;]*\));', txt, re.S))


def etype(k):
    m = re.match(r'([A-Z_0-9]+)\(', ent.get(k, ''))
    return m.group(1) if m else '?'


def refs(s):
    return re.findall(r'#(\d+)', s)


colours = {}
for k, v in ent.items():
    if v.startswith('COLOUR_RGB'):
        m = re.match(r"COLOUR_RGB\('([^']*)'", v)
        if m:
            colours[k] = m.group(1)

# MSB -> colour, by walking the style chain until a COLOUR_RGB turns up.
msb2col = {}
for k, v in ent.items():
    if not v.startswith('STYLED_ITEM'):
        continue
    r = refs(v)
    if len(r) < 2:
        continue
    psa, item = r[0], r[-1]
    if etype(item) != 'MANIFOLD_SOLID_BREP':
        continue
    seen = set()
    stack = [psa]
    col = None
    while stack:
        c = stack.pop()
        if c in seen:
            continue
        seen.add(c)
        if c in colours:
            col = colours[c]
            break
        stack.extend(refs(ent.get(c, '')))
    msb2col[item] = col
print('MSB with colour:', len(msb2col))

pts = {}
for k, v in ent.items():
    if v.startswith('CARTESIAN_POINT'):
        m = re.search(r'\(([-0-9.E+E ,]+)\)\s*\)$', v.strip())
        if m:
            try:
                c = [float(x) for x in m.group(1).split(',')]
                if len(c) == 3:
                    pts[k] = tuple(c)
            except Exception:
                pass


def bbox(root):
    """Bounding box over every CARTESIAN_POINT reachable from an entity."""
    seen = set()
    stack = [root]
    xs = []
    ys = []
    zs = []
    while stack:
        c = stack.pop()
        if c in seen:
            continue
        seen.add(c)
        if c in pts:
            p = pts[c]
            xs.append(p[0])
            ys.append(p[1])
            zs.append(p[2])
            continue
        stack.extend(refs(ent.get(c, '')))
    if not xs:
        return None
    return ((min(xs) + max(xs)) / 2, (min(ys) + max(ys)) / 2, (min(zs) + max(zs)) / 2,
            max(xs) - min(xs), max(ys) - min(ys), max(zs) - min(zs))


out = []
for k in sorted(msb2col, key=lambda k: int(k)):
    b = bbox(k)
    if b:
        out.append({'msb': int(k), 'colour': msb2col[k],
                    'cx': round(b[0], 3), 'cy': round(b[1], 3), 'cz': round(b[2], 3),
                    'dx': round(b[3], 3), 'dy': round(b[4], 3), 'dz': round(b[5], 3)})

os.makedirs(OUTDIR, exist_ok=True)
json.dump(out, open(os.path.join(OUTDIR, 'full_msb_colours.json'), 'w'), indent=1)
print('wrote output/full_msb_colours.json with', len(out), 'entries')
print('colour counts:', dict(collections.Counter(o['colour'] for o in out)))
