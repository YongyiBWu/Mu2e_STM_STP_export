"""
List every shape with all of its copies: where each sits, how it is turned, and
which part it is in the source STEP.

    "C:/Program Files/FreeCAD 1.1/bin/python.exe" scripts/list_placements.py
    "C:/Program Files/FreeCAD 1.1/bin/python.exe" scripts/list_placements.py --shape 1
    "C:/Program Files/FreeCAD 1.1/bin/python.exe" scripts/list_placements.py --csv

Answers "what is shape 8, where does it go, and what was it called in NX" in one
place. The CSVs hold the same facts but split across three files and keyed by
id, which is awkward when cross-referencing against the CAD.

Rotation
    Read from stm_placements.csv, which now records each copy's rotation as
    r11..r33 plus a readable `orientation` column -- the map from the shape's
    canonical frame to Mu2e axes, e.g. "x->x y->y z->z" for an unrotated block
    or "x->z y->y z->-x" for one turned 90 deg about y.

    This script used to re-derive that from the solid's face normals. It no
    longer does, and the difference mattered: the derivation ordered axes the
    way orthogonal_frame() returned them, while the dimensions are ordered
    canonically, so the two disagreed -- solid 0 read "x->x y->-y z->z" here
    against "x->x y->z z->-y" in the CSV. The CSV value is the one the
    round-trip verifier confirms, so it is the one reported.

    Copies with no axis-aligned rotation (wedges, 45 deg parts) carry an empty
    orientation rather than a guess.

Part id
    The `part` column of stm_placements.csv, e.g. AI-129832-A_1, taken from the
    NX document import. It is the handle for cross-referencing a row against the
    source assembly.
"""

import argparse
import collections
import csv
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
OUTDIR = os.path.join(ROOT, "output")
SIMPLE = os.path.join(ROOT, "STM_STP_files",
                      "F10269585--_1-G4 Shield House_2.stp")



def load(name):
    path = os.path.join(OUTDIR, name)
    if not os.path.exists(path):
        sys.exit("missing %s -- run scripts/extract_stm_geometry.py first"
                 % path)
    with open(path) as fh:
        return list(csv.DictReader(fh))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--shape", type=int, default=None,
                    help="only this shape_id")
    # Defaults into output/, not plots/: this is derived data keyed to the
    # CSVs, not a picture, and downstream code should read it from the same
    # place as everything else the pipeline produces.
    ap.add_argument("--csv", metavar="FILE", nargs="?", const=os.path.join(
        OUTDIR, "placements_by_shape.csv"),
                    help="also write a flat CSV of every placement "
                         "(default: output/placements_by_shape.csv)")
    # Was "--no-geometry: skip the STEP read". This script no longer reads the
    # STEP -- orientation comes from stm_placements.csv -- so that help text
    # described a saving that does not exist. All it ever gated in the end was
    # the summary tally, which is what it now honestly says.
    ap.add_argument("--no-summary", action="store_true",
                    help="omit the orientation tally printed at the end")
    args = ap.parse_args()

    shapes = {r["shape_id"]: r for r in load("stm_shapes.csv")}
    places = load("stm_placements.csv")
    bores = load("stm_bores.csv")
    nbores = collections.Counter(b["solid_id"] for b in bores)

    # Orientation comes from stm_placements.csv, not from re-deriving it here.
    #
    # This script used to recompute it from the solid's face normals in
    # orthogonal_frame() order. That is exactly the derivation the extractor
    # abandoned: it orders axes differently from the canonical frame the
    # dimensions use, so the two disagreed -- solid 0 read "x->x y->-y z->z"
    # here against "x->x y->z z->-y" in the CSV, and the CSV is the one the
    # round-trip verifier confirms. Reading the recorded value keeps this
    # listing consistent with the plots, the G4 output and the data itself.
    rot = {p["solid_id"]: (p.get("orientation") or "") for p in places}

    by_shape = collections.defaultdict(list)
    for p in places:
        by_shape[p["shape_id"]].append(p)

    wanted = ([str(args.shape)] if args.shape is not None
              else sorted(by_shape, key=int))
    if args.shape is not None and str(args.shape) not in by_shape:
        sys.exit("no shape_id %d" % args.shape)

    rows = []
    for sid in wanted:
        sh = shapes.get(sid, {})
        ps = by_shape[sid]
        size = ("%s x %s x %s mm" % (sh.get("dx"), sh.get("dy"), sh.get("dz"))
                if (sh.get("dx") or "").strip()
                else "cap %s mm^2 swept %s mm" % (sh.get("cap_area"),
                                                  sh.get("sweep_len")))
        print("=" * 78)
        print("shape %-3s %-9s %-26s %s"
              % (sid, sh.get("type", "?"), sh.get("name", "?"), size))
        print("   material %s (%s)   %d placement(s)"
              % (sh.get("material", "?"), sh.get("material_hint", "?"),
                 len(ps)))
        print()
        print("   %-6s %-18s %11s %11s %11s  %-3s %s"
              % ("solid", "part id", "x", "y", "z", "bor", "orientation"))
        for p in sorted(ps, key=lambda r: int(r["solid_id"])):
            o = rot.get(p["solid_id"], "-")
            print("   %-6s %-18s %11s %11s %11s  %-3s %s"
                  % (p["solid_id"], p.get("part", "?"), p["x"], p["y"], p["z"],
                     nbores.get(p["solid_id"], "") or "", o))
            rows.append({
                "shape_id": sid,
                "type": sh.get("type", ""),
                "name": sh.get("name", ""),
                "material": sh.get("material", ""),
                "solid_id": p["solid_id"],
                "part": p.get("part", ""),
                "x": p["x"], "y": p["y"], "z": p["z"],
                "nbores": nbores.get(p["solid_id"], 0),
                "rotY45": p.get("rotY45", ""),
                "orientation": o,
                "material_match": p.get("material_match", ""),
            })
        print()

    if args.csv:
        os.makedirs(os.path.dirname(os.path.abspath(args.csv)), exist_ok=True)
        with open(args.csv, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
            w.writeheader()
            w.writerows(rows)
        print("wrote %s (%d rows)" % (args.csv, len(rows)))

    if not args.no_summary:
        c = collections.Counter(v for v in rot.values())
        print("orientations across all %d solids:" % len(rot))
        for k, v in c.most_common():
            print("   %-28s %d" % (k, v))
    return 0


if __name__ == "__main__":
    sys.exit(main())
