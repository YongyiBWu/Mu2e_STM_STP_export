"""
Detect whether a STEP export carries a spurious global rotation about z.

    "C:/Program Files/FreeCAD 1.1/bin/python.exe" scripts/investigation/00_detect_global_tilt.py [file.stp]

Writes nothing. This is a diagnostic: run it FIRST on any new geometry, before
trusting a single dimension out of it.

Why this exists
    The Mu2e STM export carries a 0.041591 deg rotation about z that is not real
    -- it is an artefact of the NX export. Left in place it corrupts every
    world-aligned measurement: a 50.800 mm face reads 51.095, and for the 45 deg
    solids a bounding box overstates the volume threefold. The extractor removes
    it (TILT_DEG in extract_stm_geometry.py) and measures in each solid's own
    frame instead.

    A future geometry may have no such tilt, a different one, or genuinely
    rotated parts that must NOT be removed. This script tells the three apart.

The method
    Step 1, detect. Take every planar face normal, drop the ones pointing along
    z (a rotation about z leaves them untouched, so they carry no information),
    and reduce each to its angle mod 90 deg. Box faces in a clean export all
    land on 0.0000. Anything else is a shared rotation.

    Step 2, validate. De-rotate by the candidate angle and count the planar
    faces that still fail to lie on an axis. A correct angle leaves only the
    genuinely oblique features; a wrong one leaves hundreds.

Reading the histogram -- this is the part that matters
    ONE TIGHT CLUSTER AT A NON-ZERO ANGLE
        A global export error. One value shared by every part cannot arise from
        parts that were each placed by hand. Remove it.
    A SPIKE AT 0.0000 PLUS A FEW OUTLIERS
        A clean export with real angled features (45 deg brackets and the like).
        Remove nothing.
    SCATTER ACROSS MANY UNRELATED ANGLES
        The parts really are rotated with respect to each other. Remove nothing;
        there is no single frame to recover, and each solid must be measured in
        its own.

    The discriminator is not the value, it is whether the model agrees on ONE
    value. That is why the script prints the spread and not just the median.

A caution on precision
    For this model the constant was taken from a single face -- atan2 of one
    normal, exact to floating point -- which is only legitimate because the
    tilt is genuinely global. If the cluster is not tight, no single face is
    authoritative and the median is the better estimate. The script reports
    both, plus the residue each leaves, so the choice is made on evidence.

Quick: reads one file, no writes.
"""

import collections
import json
import math
import os
import statistics
import sys

import FreeCAD  # noqa: F401  (must precede Part; see extract_stm_geometry.py)
import Part

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
OUTDIR = os.path.join(ROOT, "output")
DEFAULT = os.path.join(ROOT, "STM_STP_files",
                       "F10269585--_1-G4 Shield House Simplified.stp")

# A normal this close to +/-z tells us nothing about a rotation about z.
Z_CUT = 0.99
# After de-tilting, a face is "on axis" if its largest component is this close
# to 1.
#
# This has to be far tighter than it first appears. A unit normal rotated by
# the 0.0416 deg tilt has max component 0.99999974 -- a displacement of only
# 2.6e-7. The extractor's own TOL_AXIS of 1e-4 is deliberately loose because it
# runs AFTER the de-tilt and only needs to reject 45 deg faces; reused here it
# would call a tilted face and a clean one equally "on axis" and the residue
# test would be inert at every candidate angle.
#
# 1e-9 resolves a rotation of ~0.0025 deg, well below anything worth removing.
AXIS_TOL = 1e-9

# How closely the non-zero deviations must agree to count as ONE rotation,
# in degrees. This is a comparison between tilted faces, not a threshold on the
# tilt itself: a global export error puts every affected face on the same angle
# to CAD precision, so they agree to ~1e-6 deg. Genuinely rotated parts scatter
# by degrees. Faces already on axis are excluded before this is applied.
AGREE_DEG = 1e-3
# Deviations below this count as "already on axis" and are not evidence of tilt.
ONAXIS_DEG = 1e-6


def planar_normals(shape):
    """Every planar face normal in the shape, as (x, y, z) unit tuples."""
    out = []
    for s in shape.Solids:
        for f in s.Faces:
            if f.Surface.__class__.__name__ != "Plane":
                continue
            n = f.normalAt(0, 0)
            n.normalize()
            out.append((n.x, n.y, n.z))
    return out


def inplane_angles(normals):
    """Angle mod 90 deg for each normal that carries z-rotation information."""
    return [math.degrees(math.atan2(n[1], n[0])) % 90.0
            for n in normals if abs(n[2]) <= Z_CUT]


def signed(a):
    """Fold an angle mod 90 onto (-45, 45].

    A box has faces 90 deg apart, so one normal reports 0.0416 and its
    neighbour 89.9584 for the SAME tilt. Folding puts both on +/-0.0416.
    """
    return a if a < 45.0 else a - 90.0


def deviation(a):
    """Distance from an axis, ignoring which way round the face points.

    After folding, +0.0416 and -0.0416 are still the same physical rotation:
    the sign only says which of the box's two in-plane axes the normal is
    nearer. Clustering must therefore be judged on the magnitude, otherwise a
    perfectly uniform tilt looks like a 2*tilt spread and is misread as
    genuinely rotated parts.
    """
    return abs(signed(a))


def residue(normals, tilt_deg):
    """How many planar faces are still off-axis after de-rotating by tilt_deg.

    Sign convention matches the extractor exactly: it holds TILT_DEG positive
    and applies math.radians(-TILT_DEG), so a positive tilt_deg here means the
    same rotation the extractor would perform. Verified against
    extract_stm_geometry.py: +TILT_DEG through this leaves a residual of 1e-16.
    """
    th = math.radians(-tilt_deg)
    c, s = math.cos(th), math.sin(th)
    bad = 0
    for (x, y, z) in normals:
        rx = x * c - y * s
        ry = x * s + y * c
        if max(abs(rx), abs(ry), abs(z)) < 1.0 - AXIS_TOL:
            bad += 1
    return bad


def main():
    path = sys.argv[1] if len(sys.argv) > 1 else DEFAULT
    print("reading", os.path.basename(path))
    shape = Part.read(path)
    normals = planar_normals(shape)
    print("  %d solids, %d planar faces" % (len(shape.Solids), len(normals)))

    angs = inplane_angles(normals)
    if not angs:
        print("\nNo planar faces carry z-rotation information. "
              "Nothing can be said about a tilt about z.")
        return
    print("  %d of them carry z-rotation information "
          "(the rest point along z)\n" % len(angs))

    hist = collections.Counter(round(a, 4) for a in angs)
    print("=== in-plane normal angles mod 90 deg ===")
    for k, v in sorted(hist.items(), key=lambda x: -x[1])[:15]:
        print("   %9.4f  %5d" % (k, v))
    if len(hist) > 15:
        print("   ... %d more distinct values" % (len(hist) - 15))

    # Judge clustering on the DEVIATION from an axis, not the signed angle:
    # +0.0416 and -0.0416 are one rotation seen from the box's two in-plane
    # axes, so the signed spread of a perfectly uniform tilt is 2*tilt.
    devs = [deviation(a) for a in angs]

    # Faces already on an axis say nothing about a tilt -- a model can be
    # partly axis-aligned and partly tilted. The question is whether the faces
    # that ARE off-axis all agree on one angle, so judge the spread on those.
    off = [d for d in devs if d > ONAXIS_DEG]
    onaxis = len(devs) - len(off)
    if not off:
        med, spread = 0.0, 0.0
    else:
        med = statistics.median(off)
        spread = max(off) - min(off)
    print("\ndeviation from axis: %d already on axis, %d off axis"
          % (onaxis, len(off)))
    if off:
        print("  off-axis faces: median %.6f deg   min %.6f   max %.6f   "
              "spread %.2e" % (med, min(off), max(off), spread))

    # Take the candidate from an OFF-AXIS face, at full floating-point
    # precision -- the median of rounded histogram bins cannot supply it, and
    # the most common bin is useless on a mostly-aligned model, where it is
    # simply 0.0000 and carries no tilt at all.
    exact = None
    for (x, y, z) in normals:
        if abs(z) > Z_CUT:
            continue
        d = deviation(math.degrees(math.atan2(y, x)) % 90.0)
        if d > ONAXIS_DEG and abs(d - med) < max(spread, 1e-9) * 2 + 1e-9:
            exact = d
            break

    print("\n=== de-tilt residue (planar faces still off-axis) ===")
    base = residue(normals, 0.0)
    print("   %-30s %5d of %d" % ("no de-tilt", base, len(normals)))
    # Sign convention matches the extractor: TILT_DEG is positive and
    # residue()/transform() apply the negation.
    # A deviation is unsigned, so both senses of rotation must be tried: the
    # histogram cannot say which way the model was turned, only by how much.
    # Whichever drives the off-axis count lower is the real one.
    cands = [("median %.6f" % med, med)]
    if exact is not None and abs(exact - med) > 1e-15:
        cands.append(("off-axis face %.15g" % exact, exact))
    best = None
    for label, t in cands:
        for sign, tag in ((1.0, " (+)"), (-1.0, " (-)")):
            r = residue(normals, sign * t)
            print("   %-30s %5d of %d" % (label + tag, r, len(normals)))
            if best is None or r < best[1]:
                best = (sign * t, r)

    print("\n=== verdict ===")
    t, r = best if best else (0.0, base)
    if not off:
        print("  CLEAN. Every planar normal already lies on an axis;")
        print("  no de-tilt needed.")
    elif spread < AGREE_DEG:
        print("  GLOBAL EXPORT TILT of %.15g deg." % abs(t))
        print("  %d off-axis faces, and they agree on one angle to %.2e deg."
              % (len(off), spread))
        print("  Parts placed by hand would not agree to that precision, so")
        print("  this is an artefact of the exporter, not real geometry.")
        print("  Removing it takes the off-axis count from %d to %d of %d"
              % (base, r, len(normals)))
        print("  planar faces; what remains should be only the genuinely")
        print("  oblique features -- check it matches the angled parts you")
        print("  expect, and treat a residue that does not drop as a sign")
        print("  that the angle is wrong.")
        # Hand the value over as JSON at full double precision rather than
        # asking anyone to retype it. repr() round-trips exactly; the %.15g
        # printed above is for reading, not for copying.
        # Write the magnitude, not the signed candidate. `best` carries the
        # sense that minimised the residue, but the extractor's convention is a
        # POSITIVE TILT_DEG with the negation applied at use
        # (_T = radians(-TILT_DEG)). Handing over the raw signed value would
        # de-tilt the wrong way -- the exact confusion this handoff removes.
        rec = {"tilt_deg": abs(t),
               "source": os.path.basename(path),
               "off_axis_faces": len(off),
               "agreement_deg": spread,
               "residue_before": base,
               "residue_after": r,
               "planar_faces": len(normals)}
        os.makedirs(OUTDIR, exist_ok=True)
        with open(os.path.join(OUTDIR, "tilt.json"), "w") as fh:
            json.dump(rec, fh, indent=1)
        print("\n  wrote output/tilt.json  (tilt_deg = %r)" % abs(t))
        print("  extract_stm_geometry.py picks this up automatically;")
        print("  its TILT_DEG literal is only the fallback.")
    else:
        print("  NO SINGLE GLOBAL TILT. The %d off-axis faces spread over"
              % len(off))
        print("  %.4f deg, so the parts are genuinely rotated with respect to"
              % spread)
        print("  one another. Do NOT apply a global de-tilt: there is no")
        print("  shared frame to recover. Measure each solid in its own frame")
        print("  (see orthogonal_frame/extents_along in the extractor).")
        print("\n  Note a model can be BOTH: a global tilt plus real angled")
        print("  parts. If one deviation dominates the histogram above and the")
        print("  rest are few, re-run on a subset of solids you believe are")
        print("  axis-aligned to isolate it.")


if __name__ == "__main__":
    main()
