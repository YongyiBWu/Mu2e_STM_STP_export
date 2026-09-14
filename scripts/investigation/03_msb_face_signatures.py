"""
Compute a surface-type signature per MSB entity, and check it against the kernel.

    "C:/Program Files/FreeCAD 1.1/bin/python.exe" scripts/investigation/03_msb_face_signatures.py

Writes output/full_msb_sigs.json: MSB id -> sorted (surface type, count) tuple,
e.g. [["CYLINDRICAL_SURFACE",4],["PLANE",6]] for a bored block.

Why a signature at all: centroid matching (02) is ambiguous. A signature derived
only from the STEP graph -- face count plus each face's surface type -- is
independent of coordinates, so it cannot be contaminated by the de-tilt or by
simplification moving a part slightly.

The script prints how many signature classes have identical membership counts on
the graph side and the kernel side. That agreement is what justifies the
positional bind used in 05: if the classes line up in size, the n-th member of a
class on one side is the n-th on the other.

KMAP translates FreeCAD's surface class names to STEP's vocabulary so the two
signatures are directly comparable.

Slow: reads the 74MB full assembly with the kernel.
"""

import collections
import json
import os
import re
import time

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


msb = sorted((k for k in ent if etype(k) == 'MANIFOLD_SOLID_BREP'), key=lambda k: int(k))

SURFACES = ('PLANE', 'CYLINDRICAL_SURFACE', 'CONICAL_SURFACE', 'SPHERICAL_SURFACE',
            'TOROIDAL_SURFACE', 'B_SPLINE_SURFACE', 'SURFACE_OF_LINEAR_EXTRUSION',
            'B_SPLINE_SURFACE_WITH_KNOTS')


def face_sig(k):
    """(surface type, count) for the faces of one MSB, from the graph alone."""
    shell = None
    for c in refs(ent[k]):
        if etype(c) in ('CLOSED_SHELL', 'ORIENTED_CLOSED_SHELL'):
            shell = c
            break
    if shell is None:
        return None
    types = []
    for f in refs(ent[shell]):
        if etype(f) != 'ADVANCED_FACE':
            continue
        st = '?'
        for c in refs(ent[f]):
            if etype(c) in SURFACES:
                st = etype(c)
                break
        types.append(st)
    return tuple(sorted(collections.Counter(types).items()))


sigs = {k: face_sig(k) for k in msb}
print('MSB signatures computed:', sum(1 for v in sigs.values() if v))
sc = collections.Counter(sigs.values())
print('distinct signatures:', len(sc),
      ' unique(1 solid):', sum(1 for v in sc.values() if v == 1),
      ' ambiguous groups:', sum(1 for v in sc.values() if v > 1))

t = time.time()
sh = Part.read(FULL)
print('kernel read', round(time.time() - t, 1))

KMAP = {'Plane': 'PLANE', 'Cylinder': 'CYLINDRICAL_SURFACE', 'Cone': 'CONICAL_SURFACE',
        'Sphere': 'SPHERICAL_SURFACE', 'Toroid': 'TOROIDAL_SURFACE',
        'BSplineSurface': 'B_SPLINE_SURFACE_WITH_KNOTS',
        'SurfaceOfExtrusion': 'SURFACE_OF_LINEAR_EXTRUSION'}


def ksig(s):
    c = collections.Counter(KMAP.get(f.Surface.__class__.__name__, '?') for f in s.Faces)
    return tuple(sorted(c.items()))


ks = collections.Counter(ksig(s) for s in sh.Solids)
same = sum(1 for k, v in sc.items() if ks.get(k, 0) == v)
print('signature classes matching kernel exactly:', same, 'of', len(sc))

json.dump({str(k): v for k, v in sigs.items()},
          open(os.path.join(OUTDIR, 'full_msb_sigs.json'), 'w'))
print('wrote output/full_msb_sigs.json')
