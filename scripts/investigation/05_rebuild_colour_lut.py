"""
Rebuild the dimension -> colour lookup using a validated positional bind.

    "C:/Program Files/FreeCAD 1.1/bin/python.exe" scripts/investigation/05_rebuild_colour_lut.py

Writes output/full_solid_colour.json: kernel solid index -> colour (all 224)
       output/dim_colour_lut.json:    "dx,dy,dz" (sorted, mm) -> colour

This supersedes the centroid bind in 02. The method:

    bucket MSB entities by face count   -> mb[fc]
    bucket kernel solids by face count  -> kb[fc]
    where the two buckets are the same size, the n-th of one is the n-th of the
    other, so zip them

That is only legitimate because both sides are ordered consistently (MSB by
entity id, kernel by read order) and 03 established the class sizes agree. The
script prints the class sizes and any mismatches, so the assumption is visible
rather than hidden -- classes that disagree are skipped, not guessed.

The LUT is keyed on dimensions rounded to whole mm and sorted, which makes it
orientation-free and tolerant of simplification shifting a dimension slightly.

A caution now that the tilt is known to be PIECEWISE (208 of 224 solids carry
it, 12 are square -- see stage 00): these keys come from world-aligned bounding
boxes, so they inherit whatever the solid's orientation and tilt state do to
that box. Rounding to whole mm absorbs the ~0.3mm a tilt adds to a 50.8 face,
which is why square and tilted copies of the same brick still share a key. It
does NOT absorb a genuine rotation: a brick standing on a different end gives a
different key, as shape 1 shows (143 members on one key, 3 rotated ones on two
others). That is a limitation of bbox keys, not of the tilt -- the extractor's
own dedup avoids it by measuring in each solid's own frame.
Classes that end up with more than one colour are counted as "mixed" and take
the most common; the real extractor's load_colour_lut() is the descendant of
this and marks those cases.

Slow: reads the 74MB full assembly with the kernel.
"""

import collections
import json
import os
import re

import FreeCAD  # noqa: F401
import Part

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


def shell_faces(k):
    for c in refs(ent[k]):
        if etype(c) in ('CLOSED_SHELL', 'ORIENTED_CLOSED_SHELL'):
            return [f for f in refs(ent[c]) if etype(f) == 'ADVANCED_FACE']
    return []


colours = {}
for k, v in ent.items():
    if v.startswith('COLOUR_RGB'):
        m = re.match(r"COLOUR_RGB\('([^']*)'", v)
        if m:
            colours[k] = m.group(1)

msb2col = {}
for k, v in ent.items():
    if not v.startswith('STYLED_ITEM'):
        continue
    r = refs(v)
    if len(r) < 2 or etype(r[-1]) != 'MANIFOLD_SOLID_BREP':
        continue
    seen = set()
    stack = [r[0]]
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
    msb2col[r[-1]] = col

msb = sorted((k for k in ent if etype(k) == 'MANIFOLD_SOLID_BREP'), key=lambda k: int(k))
sh = Part.read(FULL)

mb = collections.defaultdict(list)
kb = collections.defaultdict(list)
for k in msb:
    mb[len(shell_faces(k))].append(k)
for i, s in enumerate(sh.Solids):
    kb[len(s.Faces)].append(i)

print('msb classes', {k: len(v) for k, v in sorted(mb.items())})
print('ker classes', {k: len(v) for k, v in sorted(kb.items())})
bad = [fc for fc in set(list(mb) + list(kb)) if len(mb.get(fc, [])) != len(kb.get(fc, []))]
print('class size mismatches:', bad)

lut = collections.defaultdict(collections.Counter)
solidcol = {}
for fc in mb:
    if fc not in kb or len(mb[fc]) != len(kb[fc]):
        continue  # class sizes disagree: skip rather than guess
    for k, i in zip(mb[fc], kb[fc]):
        col = msb2col.get(k)
        if not col:
            continue
        solidcol[i] = col
        b = sh.Solids[i].BoundBox
        key = ','.join(str(int(round(x)))
                       for x in sorted((b.XLength, b.YLength, b.ZLength)))
        lut[key][col] += 1

print('full solids coloured:', len(solidcol), 'of', len(sh.Solids))

out = {}
mixed = 0
for k, c in lut.items():
    out[k] = c.most_common(1)[0][0]
    if len(c) > 1:
        mixed += 1

json.dump(out, open(os.path.join(OUTDIR, 'dim_colour_lut.json'), 'w'), indent=1)
json.dump({str(k): v for k, v in solidcol.items()},
          open(os.path.join(OUTDIR, 'full_solid_colour.json'), 'w'), indent=1)
print('LUT classes:', len(out), ' mixed-colour classes:', mixed)
print('colour totals:', dict(collections.Counter(solidcol.values())))
