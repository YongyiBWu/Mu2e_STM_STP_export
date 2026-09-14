"""
Inventory the simplified file: one row per solid with signature and dimensions.

    "C:/Program Files/FreeCAD 1.1/bin/python.exe" scripts/investigation/04_simplified_inventory.py

Writes output/simplified_solids.json: the 224 solids with face-type signature,
face count, centroid, bounding box and volume.

This is the side of the problem that actually gets exported. It is quick (the
simplified file is 2.2MB against the full file's 74MB) and is the input to the
cross-file match in 06.

The dimensions here are raw BOUNDING BOX values, which for this model are wrong:
the NX export carries a 0.0416 deg tilt, so a 50.800 face reads 51.095. That is
fine for the overlap matching done downstream, but it is why the real extractor
projects vertices onto recovered face normals instead. See the Dimensions note in
extract_stm_geometry.py.

The printed signature-class table is the quick view of what shapes exist.
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
                      "F10269585--_1-G4 Shield House Simplified.stp")

sh = Part.read(SIMPLE)
rows = []
for i, s in enumerate(sh.Solids):
    b = s.BoundBox
    sig = tuple(sorted(collections.Counter(
        f.Surface.__class__.__name__ for f in s.Faces).items()))
    rows.append({'i': i, 'sig': str(sig), 'nf': len(s.Faces),
                 'cx': round(b.Center.x, 3), 'cy': round(b.Center.y, 3),
                 'cz': round(b.Center.z, 3),
                 'dx': round(b.XLength, 3), 'dy': round(b.YLength, 3),
                 'dz': round(b.ZLength, 3), 'vol': round(s.Volume, 1)})

json.dump(rows, open(os.path.join(OUTDIR, 'simplified_solids.json'), 'w'), indent=1)
print('wrote output/simplified_solids.json', len(rows))

print('=== simplified signature classes ===')
for k, v in collections.Counter(r['sig'] for r in rows).most_common():
    print('   %-60s %4d' % (k, v))
