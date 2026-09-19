"""
Group the simplified solids into distinct shapes, ignoring position and rotation.

    "C:/Program Files/FreeCAD 1.1/bin/python.exe" scripts/investigation/07_dedup_groups.py

Writes output/dedup_groups.json: {key, n, ids} per shape group, largest first.

This is the prototype of the dedup that the real extractor performs, and it is
what showed the 146 identical 2x4x8 lead bricks collapse into a single
definition. The shape key is

    (face-type signature, sorted dimensions, cylinder radii, volume)

Sorting the dimensions is what makes it orientation-free: a brick laid on its
side has the same key as one standing up, so the 146 group whatever way round
they sit.

Dimensions are measured by projecting vertices onto the solid's own recovered
face normals, not read off the bounding box -- for a tilted solid a bbox is
wrong by ~0.3mm on a 50.8 face and badly wrong on the 45 deg solids. Where a
clean orthogonal frame cannot be recovered, it falls back to the bbox.

Measuring in each solid's OWN frame is also what makes this dedup immune to the
tilt being piecewise (208 of 224 solids carry it, 12 do not -- see stage 00).
A key built from world-aligned dimensions would put a square brick and a tilted
copy of the same brick in different groups. Confirmed on the real data: shape 1
has 146 members, 144 tilted and 2 square (solids 15 and 139), and they group
together. If a future export ever splits a shape that should be one, a
world-frame measurement creeping in is the first thing to suspect.

Volume is part of the key on purpose: two blocks of equal outside size but
different drilled holes hold different amounts of material and are different
parts, so they must not merge. See the "What this script does NOT do" note in
extract_stm_geometry.py.

Quick: simplified file only.
"""

import collections
import json
import os

import FreeCAD  # noqa: F401
import Part

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
OUTDIR = os.path.join(ROOT, "output")
SIMPLE = os.path.join(ROOT, "STM_STP_files",
                      "F10269585_G4_Shield_House_3.stp")

sh = Part.read(SIMPLE)


def frame(s):
    """Distinct planar face normals: the solid's own axes, if it has three."""
    ax = []
    for f in s.Faces:
        if f.Surface.__class__.__name__ != 'Plane':
            continue
        n = f.normalAt(0, 0)
        n.normalize()
        if not any(abs(abs(n.dot(u)) - 1) < 1e-6 for u in ax):
            ax.append(n)
    return ax


def shape_key(s):
    """Position/orientation-independent signature."""
    sig = tuple(sorted(collections.Counter(
        f.Surface.__class__.__name__ for f in s.Faces).items()))
    ax = frame(s)
    if len(ax) == 3 and all(abs(ax[a].dot(ax[b])) < 1e-6
                            for a in range(3) for b in range(a + 1, 3)):
        # Project every vertex onto each axis; the spread is the true extent.
        pts = [v.Point for v in s.Vertexes]
        d = [max(p.x * u.x + p.y * u.y + p.z * u.z for p in pts)
             - min(p.x * u.x + p.y * u.y + p.z * u.z for p in pts) for u in ax]
        dims = tuple(sorted(round(v, 2) for v in d))
    else:
        b = s.BoundBox
        dims = tuple(sorted(round(v, 2) for v in (b.XLength, b.YLength, b.ZLength)))
    radii = tuple(sorted({round(f.Surface.Radius, 3) for f in s.Faces
                          if f.Surface.__class__.__name__ == 'Cylinder'}))
    return (sig, dims, radii, round(s.Volume, 1))


g = collections.defaultdict(list)
for i, s in enumerate(sh.Solids):
    g[shape_key(s)].append(i)

print('=== DEDUP: %d distinct shapes from %d solids ===' % (len(g), len(sh.Solids)))
rows = sorted(g.items(), key=lambda x: -len(x[1]))
for k, v in rows[:25]:
    sig, dims, radii, vol = k
    inch = tuple(round(x / 25.4, 3) for x in dims)
    sigs = ''.join('%s%d' % (n[0], c) for n, c in sig)
    print('  x%-4d %-10s dims %-26s %-24s r=%s'
          % (len(v), sigs, str(dims), str(inch) + ' in', str(radii) if radii else '-'))

print()
print('groups with >1 instance: %d  (covering %d solids)'
      % (sum(1 for k, v in g.items() if len(v) > 1),
         sum(len(v) for v in g.values() if len(v) > 1)))
print('singletons: %d' % sum(1 for k, v in g.items() if len(v) == 1))

json.dump([{'key': str(k), 'n': len(v), 'ids': v} for k, v in rows],
          open(os.path.join(OUTDIR, 'dedup_groups.json'), 'w'), indent=1)
print('wrote output/dedup_groups.json')
