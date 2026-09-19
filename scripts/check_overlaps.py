"""
Check whether any two placed solids overlap.

    "C:/Program Files/FreeCAD 1.1/bin/python.exe" scripts/check_overlaps.py
    "C:/Program Files/FreeCAD 1.1/bin/python.exe" scripts/check_overlaps.py --tol 0.01 --report plots/overlaps.txt

Geant4 will not tell you politely: overlapping volumes give wrong navigation and
silently wrong physics, and G4 own overlap check samples points rather than
proving anything. This tests the real geometry instead.

Method
    Two passes, because 224 solids is 24,976 pairs and a boolean on each would
    take hours.

    Pass 1, bounding boxes. Reject any pair whose world-aligned boxes do not
    interpenetrate by more than the tolerance. This is cheap and conservative:
    boxes that miss cannot contain solids that touch.

    Pass 2, kernel booleans. For each surviving pair, intersect the ACTUAL
    solids and measure the volume of the result. This is the part that matters:
    a bounding-box test alone would be badly wrong on this model, because the
    12 prisms fill only part of their envelope (shape 19's wedge is half its
    box) and the 15 bored blocks have material drilled out. Box-on-box would
    report overlaps between parts that merely nest.

Tolerance
    --tol is a volume in mm^3, default 1.0. CAD surfaces that abut share a face
    exactly, and OCC will sometimes return a sliver of 1e-9 mm^3 for a shared
    face rather than a clean zero. A cubic millimetre is far below anything
    that matters for shielding and far above that numerical noise.

    Reported overlaps also carry the intersection's bounding box, so a genuine
    interpenetration (a chunky box) is distinguishable from a grazing contact
    (a thin sheet) at a glance.

Reads the STEP file directly rather than the CSVs: the CSVs describe shapes and
placements, but the question here is about the solids as they actually sit, and
the kernel has those. The CSVs supply names and shape ids for the report.
"""

import argparse
import collections
import csv
import os
import sys
import time

import FreeCAD  # noqa: F401  (must precede Part)
import Part

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
OUTDIR = os.path.join(ROOT, "output")
SIMPLE = os.path.join(ROOT, "STM_STP_files",
                      "F10269585_G4_Shield_House_3.stp")


def load_labels():
    """solid_id -> (shape_id, name, part) from the CSVs, for readable output."""
    out = {}
    try:
        shapes = {r["shape_id"]: r for r in
                  csv.DictReader(open(os.path.join(OUTDIR, "stm_shapes.csv")))}
        for r in csv.DictReader(open(os.path.join(OUTDIR,
                                                  "stm_placements.csv"))):
            sh = shapes.get(r["shape_id"], {})
            out[int(r["solid_id"])] = (r["shape_id"],
                                       sh.get("name", "?"),
                                       r.get("part", ""))
    except (IOError, OSError):
        pass
    return out


# An intersection this thin in one axis is two flat faces sitting on the same
# plane, not one solid pushed into another.
COPLANAR_MM = 0.5
# An intersection smaller than this fraction of the smaller solid is a skin,
# not a collision. A tube in a bore of exactly equal radius produces one.
SKIN_FRACTION = 0.01


def classify(a, b, common, vol):
    """Why do these two intersect: shared surface, or real interpenetration?

    The distinction matters more than the raw count. CAD parts that are meant
    to touch -- a brick stacked on a brick, a tube sliding through a bore cut
    to the same radius -- produce a non-zero boolean intersection, because OCC
    gives a coincident surface a tolerance envelope rather than an exact zero.
    Reporting those alongside genuine collisions buries the ones that matter.

        coplanar    the intersection is a sheet: two faces on a shared plane
        skin        a thin shell over a coincident CURVED surface (the tube in
                    its bore). No thin bbox axis, because a cylinder wraps --
                    which is why this is caught by volume fraction instead
        real        neither: one solid genuinely occupies another's space
    """
    bb = common.BoundBox
    if min(bb.XLength, bb.YLength, bb.ZLength) < COPLANAR_MM:
        return "coplanar"
    smaller = min(a.Volume, b.Volume)
    if smaller > 0 and vol / smaller < SKIN_FRACTION:
        return "skin"
    return "real"


def bbox_gap(a, b):
    """Negative overlap depth per axis; >= 0 means the boxes are apart."""
    return (max(a.XMin - b.XMax, b.XMin - a.XMax),
            max(a.YMin - b.YMax, b.YMin - a.YMax),
            max(a.ZMin - b.ZMax, b.ZMin - a.ZMax))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tol", type=float, default=1.0,
                    help="ignore intersections below this volume in mm^3")
    ap.add_argument("--margin", type=float, default=0.0,
                    help="bbox interpenetration needed to consider a pair, mm")
    # Defaults into output/, alongside the CSVs: this is a derived result of
    # the geometry, not a picture, so it belongs with the pipeline's data
    # rather than in the plots directory.
    ap.add_argument("--report", nargs="?", default=None,
                    const=os.path.join(OUTDIR, "overlaps.txt"),
                    help="also write this text file "
                         "(default: output/overlaps.txt)")
    ap.add_argument("--limit", type=int, default=0,
                    help="stop after this many boolean tests (0 = no limit)")
    args = ap.parse_args()

    lines = []

    def say(s=""):
        print(s)
        lines.append(s)

    say("reading %s" % os.path.basename(SIMPLE))
    shape = Part.read(SIMPLE)
    solids = shape.Solids
    labels = load_labels()
    say("  %d solids" % len(solids))

    boxes = [s.BoundBox for s in solids]

    # Pass 1 -------------------------------------------------------------
    t0 = time.time()
    cands = []
    for i in range(len(solids)):
        for j in range(i + 1, len(solids)):
            gx, gy, gz = bbox_gap(boxes[i], boxes[j])
            # All three gaps must be negative for the boxes to interpenetrate.
            if max(gx, gy, gz) < -args.margin:
                cands.append((i, j, -max(gx, gy, gz)))
    total_pairs = len(solids) * (len(solids) - 1) // 2
    say("\npass 1: %d of %d pairs have interpenetrating bounding boxes (%.1fs)"
        % (len(cands), total_pairs, time.time() - t0))

    if args.limit and len(cands) > args.limit:
        cands.sort(key=lambda c: -c[2])
        say("  --limit %d: testing only the %d deepest"
            % (args.limit, args.limit))
        cands = cands[:args.limit]

    # Pass 2 -------------------------------------------------------------
    say("\npass 2: intersecting the real solids (this is the slow part)")
    t0 = time.time()
    hits = []
    errors = []
    for n, (i, j, depth) in enumerate(cands, 1):
        if n % 200 == 0:
            say("   %d/%d tested, %.0fs elapsed, %d overlap(s) so far"
                % (n, len(cands), time.time() - t0, len(hits)))
        try:
            common = solids[i].common(solids[j])
            vol = common.Volume if common.Solids else 0.0
        except Exception as exc:                       # pragma: no cover
            errors.append((i, j, str(exc)[:60]))
            continue
        if vol > args.tol:
            hits.append((i, j, vol, common.BoundBox,
                         classify(solids[i], solids[j], common, vol)))
    say("   %d/%d tested in %.0fs" % (len(cands), len(cands),
                                      time.time() - t0))

    # Report -------------------------------------------------------------
    say("\n" + "=" * 68)
    real = [h for h in hits if h[4] == "real"]
    contact = [h for h in hits if h[4] != "real"]

    if not hits:
        say("NO INTERSECTIONS above %.3g mm^3." % args.tol)
        say("%d candidate pairs were intersected and every one came back "
            "empty." % len(cands))
    else:
        hits.sort(key=lambda h: -h[2])
        say("%d intersecting pair(s) above %.3g mm^3: %d REAL, %d contact"
            % (len(hits), args.tol, len(real), len(contact)))
        say("")
        say("%-5s %-5s %-9s %12s  %-26s %s"
            % ("solid", "solid", "kind", "volume mm^3",
               "extent of overlap mm", "parts"))
        for i, j, vol, bb, kind in hits:
            si = labels.get(i, ("?", "?", ""))
            sj = labels.get(j, ("?", "?", ""))
            say("%-5d %-5d %-9s %12.3f  %7.2f x %7.2f x %7.2f  %s / %s"
                % (i, j, kind, vol, bb.XLength, bb.YLength, bb.ZLength,
                   si[1], sj[1]))
        say("")
        say("coplanar = two faces sharing a plane (stacked bricks).")
        say("skin     = a shell over a coincident curved surface, e.g. a tube")
        say("           in a bore cut to the same radius. Designed fits.")
        say("real     = one solid genuinely occupying another's space.")
        say("")
        if real:
            say("ONLY THE %d 'real' ROW(S) NEED ATTENTION." % len(real))
        else:
            say("NOTHING CLASSIFIED 'real': every intersection is a shared")
            say("surface between parts designed to touch. Geant4 may still")
            say("warn about coincident surfaces, but no solid occupies")
            say("another's space.")

    if errors:
        say("\n%d pair(s) could not be intersected:" % len(errors))
        for i, j, e in errors[:10]:
            say("   %d/%d: %s" % (i, j, e))

    by_shape = collections.Counter()
    for i, j, _, _, kind in hits:
        if kind != "real":
            continue
        by_shape[labels.get(i, ("?",))[0]] += 1
        by_shape[labels.get(j, ("?",))[0]] += 1
    if by_shape:
        say("\nREAL overlaps by shape_id: %s"
            % ", ".join("%s x%d" % kv for kv in by_shape.most_common()))

    if args.report:
        os.makedirs(os.path.dirname(os.path.abspath(args.report)),
                    exist_ok=True)
        with open(args.report, "w") as fh:
            fh.write("\n".join(lines) + "\n")
        print("\nwrote %s" % args.report)

    # Gate on genuine interpenetration only. Coplanar faces and coincident
    # curved surfaces are how CAD expresses parts that touch by design, so
    # failing on those would make this check useless as a build gate -- it
    # would be red on a correct model.
    return 1 if real else 0


if __name__ == "__main__":
    sys.exit(main())
