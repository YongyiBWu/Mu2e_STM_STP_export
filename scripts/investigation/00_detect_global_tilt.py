"""
Detect whether a STEP export carries a spurious rotation about z, and whether
it applies to the whole model or only to some solids.

    "C:/Program Files/FreeCAD 1.1/bin/python.exe" scripts/investigation/00_detect_global_tilt.py [file.stp]

Writes output/tilt.json when -- and only when -- a single global angle is
actually correct for the model. Run this FIRST on any new geometry, before
trusting a single dimension out of it.

Why this exists
    The Mu2e STM export carries a 0.041591 deg rotation about z that is not
    real -- an artefact of the NX export. Left in place it corrupts every
    world-aligned measurement: a 50.800 mm face reads 51.095, and for the 45
    deg solids a bounding box overstates the volume threefold.

    But the tilt is NOT global. In this model 208 of 224 solids carry it and 12
    are already square in the CAD -- among them the SSC cluster on the beamline.
    De-rotating all 224 fixes the 208 and rotates the 12 off their axes,
    introducing the very error it removes elsewhere.

    So the first question is not "what is the angle" but "which solids have
    it". This script answers that before it answers anything else.

The method, in order
    Step 1, candidate angle. Take every planar face normal, drop the ones
    pointing along z (a rotation about z leaves them untouched, so they carry
    no information), and reduce each to its deviation from the nearest axis
    mod 90 deg. The off-axis faces' shared value is the candidate.

    Step 2, PER SOLID. Classify each solid on its own normals: square, tilted,
    oblique or indeterminate. This is the step a pooled histogram cannot do,
    and it is why the earlier version of this script missed the split: "every
    solid tilted by X" and "most solids tilted by X, some square" produce the
    same single tight cluster when the faces are averaged together.

    Step 3, validate. De-rotate by the candidate and count the planar faces
    still off-axis. A correct angle leaves only the genuinely oblique features.

Reading the verdict
    PIECEWISE TILT
        Some solids tilted, some square. A single rotation is wrong: de-tilt
        per solid, leaving the square ones alone. No tilt.json is written,
        because handing one global angle to the extractor is what causes the
        damage. Cross-check the split against the full assembly first -- if
        both files agree solid-for-solid it is real CAD structure.
    GLOBAL EXPORT TILT
        Every judgeable solid carries the same angle. One rotation is correct;
        tilt.json is written for the extractor.
    CLEAN
        Everything already on axis. Remove nothing.
    NO SINGLE TILT
        Off-axis angles scatter over degrees: the parts are genuinely rotated
        with respect to one another and there is no shared frame to recover.

    The discriminator is never the value alone. It is whether the solids agree
    -- which is why the per-solid table is printed above the verdict.

Why the damage is easy to miss
    2x the tilt displaces a unit normal by ~1e-6, far inside a 1e-4 axis-
    alignment test, and the box paths use bounding-box dimensions rather than
    the direction vectors. A wrong de-tilt therefore shows up only where a
    transformed direction reaches an output file -- bore axes and prism basis
    vectors -- long after the positions have been silently shifted.

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
                       "F10269585--_1-G4 Shield House_2.stp")

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


def solid_state(solid, tilt_deg):
    """How this ONE solid sits: 'square', 'tilted', 'oblique' or 'z-parallel'.

    Judged on the solid's own in-plane normals, independently of every other
    solid. That independence is the whole point: a pooled histogram averages
    the populations together and reports a single global tilt even when the
    model does not have one.

        square      every in-plane normal already on an axis
        tilted      every in-plane normal at tilt_deg from an axis
        oblique     carries some other angle (genuine 45 deg features)
        z-parallel  every direction on the solid is parallel to z, so a
                    rotation ABOUT z maps each of them to itself and the
                    solid's own geometry cannot say whether it was rotated

    'z-parallel' is a statement about what is observable, not a measurement.
    The two tubes are the case: both cylinder axes and both cap normals have
    an in-plane component of 6e-17 or exactly 0. Classifying them on their
    sweep would return 'square' every time as an artefact of that axis being
    the rotation axis -- a confident-looking non-answer. mating_bore_state()
    below infers their real state from what they are seated in instead.
    """
    devs = set()
    for f in solid.Faces:
        if f.Surface.__class__.__name__ != "Plane":
            continue
        n = f.normalAt(0, 0)
        n.normalize()
        if abs(n.z) > Z_CUT:
            continue
        devs.add(deviation(math.degrees(math.atan2(n.y, n.x)) % 90.0))
    if not devs:
        return "z-parallel"
    on = [d for d in devs if d <= ONAXIS_DEG]
    at = [d for d in devs if abs(d - tilt_deg) <= AGREE_DEG]
    other = [d for d in devs if d > ONAXIS_DEG and abs(d - tilt_deg) > AGREE_DEG]
    if other:
        return "oblique"
    if on and not at:
        return "square"
    if at and not on:
        return "tilted"
    return "mixed"


# How close two axes must be, in mm, to count as one seated in the other.
# Loose on purpose: most neighbours here sit 0.5-0.93 mm off, so this is
# "collinear within a millimetre", not a precision fit.
MATE_PERP_TOL = 2.0


def mating_bore_state(solids, index, states):
    """State of a z-parallel solid, inferred from the bores it sits in.

    A z-parallel solid cannot be judged on its own faces, but it is not
    floating free: both tubes here run through bores in the surrounding walls,
    and a tube seated in a tilted wall's bore is tilted with it.

    So: find every cylindrical face elsewhere that is collinear with this
    solid's own axis, and report the state of those solids if they agree.
    Returns (state, note) or (None, reason) when nothing can be concluded.

    This is an INFERENCE from mating geometry, never a measurement of the
    solid, which is why the caller records it with by="mating-bore". The
    offsets are also loose -- most neighbours here sit 0.5-0.93 mm off, so
    "collinear within MATE_PERP_TOL" is the honest reading, not a precision
    fit. Only solid 222 seats tube 157 at 0.00000 mm.
    """
    me = solids[index]
    cyl = [f for f in me.Faces if f.Surface.__class__.__name__ == "Cylinder"]
    if not cyl:
        return None, "no cylindrical face to match on"
    axis = cyl[0].Surface.Axis
    origin = cyl[0].Surface.Center
    bb = me.BoundBox

    seen = {}
    for j, s in enumerate(solids):
        if j == index:
            continue
        b = s.BoundBox
        if (b.XMin > bb.XMax + 5 or b.XMax < bb.XMin - 5
                or b.YMin > bb.YMax + 5 or b.YMax < bb.YMin - 5
                or b.ZMin > bb.ZMax + 5 or b.ZMax < bb.ZMin - 5):
            continue
        for f in s.Faces:
            if f.Surface.__class__.__name__ != "Cylinder":
                continue
            fa = f.Surface.Axis
            if abs(abs(fa.dot(axis)) - 1.0) > 1e-6:
                continue
            fc = f.Surface.Center
            d = FreeCAD.Vector(fc.x - origin.x, fc.y - origin.y,
                               fc.z - origin.z)
            if d.sub(axis * d.dot(axis)).Length <= MATE_PERP_TOL:
                seen[j] = states.get(j)
                break
    if not seen:
        return None, "not seated in any bore"

    # 'mixed' counts as carrying the tilt: such a solid straddles the split.
    kinds = {st for st in seen.values() if st is not None}
    if kinds <= {"tilted", "mixed"} and kinds:
        return "tilted", "seated in %d bore(s), all tilted" % len(seen)
    if kinds == {"square"}:
        return "square", "seated in %d bore(s), all square" % len(seen)
    return None, "seated in bores of mixed state: %s" % ", ".join(sorted(kinds))


def per_solid_report(shape, tilt_deg):
    """Classify every solid and print the table. Returns {state: [indices]}.

    This runs BEFORE any global verdict, because whether a single angle is
    even meaningful depends on the answer. A model where every solid is
    'tilted' has a global tilt; one that splits into 'tilted' and 'square'
    has a per-solid tilt, and de-rotating all of it corrupts the square parts.
    """
    solids = shape.Solids
    groups = collections.defaultdict(list)
    rows = []

    # Pass 1: every solid on its own faces.
    own = {}
    for i, s in enumerate(solids):
        own[i] = solid_state(s, tilt_deg)

    # Pass 2: a z-parallel solid has no usable faces, so ask what it is seated
    # in. Runs second because it reads pass 1's verdicts for the neighbours.
    inferred = {}
    for i, st in own.items():
        if st != "z-parallel":
            continue
        got, why = mating_bore_state(solids, i, own)
        if got is not None:
            inferred[i] = (got, why)

    for i, s in enumerate(solids):
        st = own[i]
        by = "face-normals"
        note = None
        if i in inferred:
            st, note = inferred[i]
            by = "mating-bore"
        elif st == "z-parallel":
            note = "every direction parallel to z; rotation about z unobservable"
        groups[st].append(i)
        b = s.BoundBox
        # x/y/z, not cx/cy/cz: scripts/plot_tilt_state.py pairs each centre
        # with its extent as r[k] and r["d" + k], so the names have to match.
        rows.append({"i": i, "state": st,
                     # How that state was reached: measured from the solid's
                     # own face normals, or inferred from the bores it sits
                     # in. A consumer that treats an inference as a
                     # measurement would be wrong to, so say which.
                     "by": by, "note": note,
                     "x": round(b.Center.x, 3), "y": round(b.Center.y, 3),
                     "z": round(b.Center.z, 3),
                     "dx": round(b.XLength, 3), "dy": round(b.YLength, 3),
                     "dz": round(b.ZLength, 3),
                     "vol": round(s.Volume, 1)})

    # Write the per-solid table, not just print it. A per-solid de-tilt needs to
    # know WHICH solids to rotate, so this classification is data the rest of
    # the chain consumes -- scripts/plot_tilt_state.py draws from it -- rather
    # than something to re-derive by hand each time. Centres are CAD-frame here
    # on purpose: this file describes the INPUT, before any transform.
    os.makedirs(OUTDIR, exist_ok=True)
    with open(os.path.join(OUTDIR, "tilt_state.json"), "w") as fh:
        json.dump(rows, fh, indent=1)

    total = len(shape.Solids)
    print("=== per-solid tilt state (candidate %.6f deg) ===" % tilt_deg)
    for state in ("tilted", "square", "oblique", "mixed", "z-parallel"):
        ids = groups.get(state)
        if not ids:
            continue
        shown = ", ".join(str(i) for i in ids[:10])
        if len(ids) > 10:
            shown += ", ... (%d more)" % (len(ids) - 10)
        print("  %-11s %4d of %d   %s" % (state, len(ids), total, shown))
    if inferred:
        print("  (%d solid(s) classified from mating bores, not from their"
              " own faces:)" % len(inferred))
        for i in sorted(inferred):
            st, why = inferred[i]
            print("     solid %-4d -> %-8s  %s" % (i, st, why))
    return groups


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

    t, r = best if best else (0.0, base)

    # Classify each solid on its own BEFORE pronouncing on the model as a
    # whole. A pooled histogram cannot tell "every solid is tilted by X" from
    # "most solids are tilted by X and some are already square" -- both give
    # one tight off-axis cluster. The difference decides whether a single
    # global de-tilt is correct or actively destructive.
    print()
    # Always classify, even when nothing is off-axis. A clean model still owes
    # the chain a tilt_state.json saying so -- "every solid square" is a
    # result, not an absence of one, and plot_tilt_state.py needs the file
    # either way. With no off-axis faces the candidate is 0.0, which makes
    # every solid classify as 'square'.
    groups = per_solid_report(shape, abs(t))
    n_tilted = len(groups.get("tilted", []))
    n_square = len(groups.get("square", []))
    piecewise = bool(off) and n_tilted and n_square

    print("\n=== verdict ===")
    if not off:
        print("  CLEAN. Every planar normal already lies on an axis;")
        print("  no de-tilt needed.")
    elif piecewise:
        print("  PIECEWISE TILT of %.15g deg -- NOT global." % abs(t))
        print("  %d solids carry it, %d are already square in the CAD."
              % (n_tilted, n_square))
        print()
        print("  A single de-tilt is WRONG for this model. It would correct")
        print("  the %d tilted solids and rotate the %d square ones off their"
              % (n_tilted, n_square))
        print("  axes, introducing exactly the error it removes elsewhere.")
        print("  The damage is easy to miss: %.15g deg displaces a unit normal"
              % abs(t))
        print("  by only %.2e, well inside a 1e-4 axis test, and box paths use"
              % (1.0 - math.cos(math.radians(abs(t)))))
        print("  bounding-box dimensions rather than the direction vectors.")
        print("  It surfaces only where a transformed direction reaches an")
        print("  output -- bore axes and prism basis vectors.")
        print()
        print("  De-tilt PER SOLID: apply the rotation to the tilted set only.")
        print("  The square solids listed above must be left alone.")
        print()
        print("  Cross-check the split against the full assembly before")
        print("  acting: if both files agree solid-for-solid, the split is")
        print("  real CAD structure and not an artefact of simplification.")
    elif spread < AGREE_DEG:
        print("  GLOBAL EXPORT TILT of %.15g deg." % abs(t))
        print("  %d off-axis faces, and they agree on one angle to %.2e deg."
              % (len(off), spread))
        print("  Every solid that can be judged carries it, so one rotation")
        print("  is correct for the whole model.")
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
