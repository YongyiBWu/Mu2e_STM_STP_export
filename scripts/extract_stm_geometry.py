"""
Extract the STM shield-house geometry from the NX STEP export into a
deduplicated, classified inventory in Mu2e coordinates.

Run with FreeCAD's bundled interpreter (FreeCAD must be imported before Part):

    "C:/Program Files/FreeCAD 1.1/bin/python.exe" scripts/extract_stm_geometry.py

Inputs
    STM_STP_files/F10269585--_1-G4 Shield House Simplified.stp   (224 solids)
    STM_STP_files/F10258491--_1-Shield House Square.stp          (390 solids)

Outputs (output/)
    stm_shapes.csv      one row per distinct shape: the block definitions
    stm_placements.csv  one row per solid: which shape goes where
    stm_bores.csv       cylindrical cuts, for the G4SubtractionSolid cases

Coordinate transform (CAD -> Mu2e)
    1. De-tilt by +0.041591 deg about z. The NX export carries a spurious
       rotation: planar normals cluster at 0.0416 and 89.9584 deg mod 90, one
       global value, not per-part scatter. It is a CAD error, so it is removed.
       Afterwards only 28 of 1364 planar faces are off-axis, and those are
       genuine 45 deg features.
    2. CAD x -> -x, y -> y, z -> -z, i.e. 180 deg about y.
    3. Subtract the anchor so it lands on _STMShieldingRef.

    Anchor: face 21 of AI-129832 (the tungsten spot-size collimator), its
    low-z face, CAD (-354.292, -23.999, 370.195). Offline places the SSC at
    _STMShieldingRef + (0,0,Wdepth_f/2), so the SSC front face sits exactly on
    _STMShieldingRef. Subtracting the anchor in all three coordinates also puts
    the two SSC bores at x = +/-40.64, reproducing Offline's
    stm.STM_SSC.offset_Spot = 40.64 -- an independent check that the fit is right.

Dimensions
    Box sizes come from projecting vertices onto the recovered face normals,
    never from the bounding box. With the model tilted a bbox reads 51.095 for
    a face that is really 50.800, and for the 45 deg solids it is wrong by a
    factor of three.

Materials
    The simplified STEP carries no material data: no material entities at all,
    and a single colour for the whole file. The full assembly is colour-coded,
    so material is recovered from there in two hops:

        simplified solid --(bounding-box overlap)--> full solid
        full solid       --(face-count bucket)-----> STEP colour

    The second hop is the awkward one and is explained at _step_solid_colours.
    The value written to the CSV is the CAD colour name, e.g. "granite gray",
    as a placeholder to be replaced by a real Geant4 material; COLOUR_HINT
    carries a guess alongside it. Every solid also gets a material_match score
    so a material inferred from a loose match is visible rather than implied.

Names
    Suggestions only, meant to be edited. Blocks that repeat and are scattered
    through the house are named for what they are (LeadBrick2x4x8PV); panels
    that stay in one place follow the existing Offline scheme,
    {Material}{Side}wall{N}PV (CopperLwallPV, LeadTwall1PV, BPRwall2PV). The
    side letter comes from a position heuristic and is not authoritative.

What this script does NOT do
    It does not merge shapes that differ only by rotation *and* by drilled
    holes. Sorted dimensions already make the shape key orientation-free, so
    the 146 bricks collapse into one group whatever way round they sit. Blocks
    of equal outside size but different bores stay separate on purpose: they
    hold different amounts of material (e.g. 1048772 vs 945809 mm^3) and are
    different parts.
"""

import collections
import csv
import json
import math
import os
import re

# FreeCAD must be imported before Part. Importing Part first fails with
# "No module named 'Part'": it is FreeCAD that puts the kernel DLLs on the
# search path. This is also why the script has to run under FreeCAD's own
# interpreter rather than a system Python.
import FreeCAD
import Part

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
OUTDIR = os.path.join(ROOT, "output")

SIMPLE = os.path.join(ROOT, "STM_STP_files",
                      "F10269585--_1-G4 Shield House Simplified.stp")
FULL = os.path.join(ROOT, "STM_STP_files",
                    "F10258491--_1-Shield House Square.stp")

def _load_tilt():
    """The de-tilt angle, preferring the value measured from the geometry.

    scripts/investigation/00_detect_global_tilt.py measures the tilt and writes
    output/tilt.json. Reading it back keeps full double precision and removes
    the transcription step: on a new export, run the detector and this picks the
    new value up. The literal below is the fallback, and is the value measured
    from F10269585 (it agrees with the detector to 3e-14 deg).

    A tilt.json measured from a DIFFERENT file than the one being extracted is
    ignored with a warning -- a tilt is a property of one export, not a setting.
    """
    fallback = 0.041591442757606956
    path = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        os.pardir, "output", "tilt.json")
    try:
        with open(path) as fh:
            rec = json.load(fh)
    except (IOError, OSError, ValueError):
        return fallback
    src = rec.get("source")
    if src and src != os.path.basename(SIMPLE):
        print("warning: output/tilt.json was measured from %s, not %s"
              % (src, os.path.basename(SIMPLE)))
        print("         ignoring it and using the built-in value")
        return fallback
    tilt = rec.get("tilt_deg")
    if not isinstance(tilt, float):
        return fallback
    return tilt


TILT_DEG = _load_tilt()             # removed, see module docstring

# Anchor in the de-tilted, y-flipped frame. With this value face 21 of
# AI-129832 lands at (0.0002, -0.0002, 0.0003) and the SSC bores at
# x = +/-40.6400, reproducing Offline's stm.STM_SSC.offset_Spot = 40.64.
ANCHOR = (354.3093, -23.7418, -370.1950)

# Comparing unit directions: two face normals are the same axis if |n.u| = 1.
TOL_UNIT = 1e-6
# Deciding "does this normal point along x, y or z". Loose enough to survive
# the de-tilt leaving ~1e-5 residue, tight enough to reject a 45 deg face.
TOL_AXIS = 1e-4
# Relative volume agreement, used to confirm a solid really is the primitive
# its faces suggest. A box whose volume matches dx*dy*dz to 1e-4 is a box.
TOL_VOL = 1e-4

# Colour name -> the Geant4 material it most likely stands for. Deliberately
# left as a hint: edit the material column rather than trusting these.
COLOUR_HINT = {
    "granite gray": "G4_Pb?",
    "medium steel": "RackSteel?",
    "medium maroon": "G4_Cu?",
    "strong ice": "BP?",
    "silver gray": "G4_Al?",
    "medium ice": "Polyethylene?",
    "pale sky": "?",
    "powder gray": "G4_W?",
    "ash gray": "?",
    "pale red": "?",
}

MATERIAL_WORD = {
    "granite gray": "Lead",
    "medium steel": "Steel",
    "medium maroon": "Copper",
    "strong ice": "BP",
    "silver gray": "Aluminum",
    "medium ice": "Poly",
    "powder gray": "Tungsten",
}


# ---------------------------------------------------------------- transform

# Rotating by -TILT_DEG undoes the tilt the exporter baked in.
_T = math.radians(-TILT_DEG)


def transform(x, y, z):
    """CAD point -> Mu2e point relative to _STMShieldingRef.

    Three steps in one expression: de-tilt about z, flip 180 deg about y
    (which is what negating x and z does), then shift onto the anchor.
    """
    xr = x * math.cos(_T) - y * math.sin(_T)
    yr = x * math.sin(_T) + y * math.cos(_T)
    return (-xr - ANCHOR[0], yr - ANCHOR[1], -z - ANCHOR[2])


def transform_dir(x, y, z):
    """Same rotation as transform(), without the translation.

    Directions must not be shifted, so normals and axes come through here.
    """
    xr = x * math.cos(_T) - y * math.sin(_T)
    yr = x * math.sin(_T) + y * math.cos(_T)
    return (-xr, yr, -z)


def is_axis_aligned(n):
    """Does this direction lie along x, y or z?

    One component near +/-1 means the other two are near zero, since n is a
    unit vector. Used to tell a de-tilted face from a genuine 45 deg one.
    """
    return max(abs(n[0]), abs(n[1]), abs(n[2])) > 1.0 - TOL_AXIS


# ------------------------------------------------------------- geometry aid

def unique_normals(solid):
    """The distinct planar face directions of a solid.

    Opposite faces of a box share an axis, so n and -n count once: a box
    returns three directions, not six. The count is a cheap shape fingerprint
    (3 -> box-like, 4+ -> wedge or skewed prism).
    """
    axes = []
    for f in solid.Faces:
        if f.Surface.__class__.__name__ != "Plane":
            continue
        n = f.normalAt(0, 0)
        n.normalize()
        if not any(abs(abs(n.dot(u)) - 1.0) < TOL_UNIT for u in axes):
            axes.append(n)
    return axes


def orthogonal_frame(solid):
    """The solid's own three axes, if it has three mutually square ones.

    This is the local coordinate system a rotated box is built in. Returning
    it lets dimensions be measured in the solid's frame instead of the world's.
    None means the solid is not box-like.
    """
    axes = unique_normals(solid)
    if len(axes) != 3:
        return None
    for a in range(3):
        for b in range(a + 1, 3):
            if abs(axes[a].dot(axes[b])) > TOL_UNIT:
                return None
    return axes


def extents_along(solid, axes):
    """Size and mid-point along each of the solid's own axes.

    Projecting every vertex onto an axis and taking the span gives the true
    edge length. This is the reason the script does not use BoundBox: with the
    model tilted, the world-aligned box of a 50.800 face reads 51.095, and for
    a 45 deg solid it overstates the volume threefold.
    """
    pts = [v.Point for v in solid.Vertexes]
    dims, mids = [], []
    for u in axes:
        pr = [p.x * u.x + p.y * u.y + p.z * u.z for p in pts]
        dims.append(max(pr) - min(pr))
        mids.append((max(pr) + min(pr)) / 2.0)
    return dims, mids


def cylinders(solid):
    """Cylindrical faces: bores if the solid is a block, walls if it is a tube."""
    return [f for f in solid.Faces
            if f.Surface.__class__.__name__ == "Cylinder"]


def oblique(solid):
    """True if the solid keeps off-axis faces after the de-tilt (45 deg parts)."""
    for f in solid.Faces:
        if f.Surface.__class__.__name__ != "Plane":
            continue
        if not is_axis_aligned(transform_dir(*f.normalAt(0, 0))):
            return True
    return False


# ------------------------------------------------------------ classification

def classify(solid):
    """(kind, params) with lengths in mm and positions in Mu2e coordinates.

    Tried in order, most specific first. Each test is confirmed against the
    solid's volume rather than trusting face counts alone, so a block that
    merely looks like a box does not silently become a G4Box:

        BOX       6 planar faces, 3 square axes, volume == dx*dy*dz
        TUBE      2 coaxial cylinders + 2 caps, volume == annulus
        BOX_HOLE  planar block with cylindrical bores -> G4SubtractionSolid
        PRISM     all-planar but not a box: wedges, chamfered blocks
        OTHER     anything left, reported so nothing is dropped in silence

    For this model the tally is 195 / 2 / 15 / 12 / 0.
    """
    surf = {f.Surface.__class__.__name__ for f in solid.Faces}
    cyls = cylinders(solid)

    # --- BOX: six planes meeting at right angles ------------------------
    if surf == {"Plane"} and len(solid.Faces) == 6:
        axes = orthogonal_frame(solid)
        if axes is not None:
            dims, mids = extents_along(solid, axes)
            product = dims[0] * dims[1] * dims[2]
            # Volume check: three pairs of parallel planes could still be a
            # sheared box, whose volume would fall short of the product.
            if product > 0 and abs(solid.Volume / product - 1.0) < TOL_VOL:
                # Centre of the solid, rebuilt from the mid-point along each
                # of its own axes.
                c = FreeCAD.Vector(0, 0, 0)
                for u, m in zip(axes, mids):
                    c = c + u * m
                # Re-label the three measured lengths as dx/dy/dz by asking
                # which world axis each of the solid's axes points along once
                # transformed. Without this a rotated box would report its
                # height as a width.
                sx = sy = sz = None
                for u, d in zip(axes, dims):
                    n = transform_dir(u.x, u.y, u.z)
                    if is_axis_aligned(n):
                        k = max(range(3), key=lambda i: abs(n[i]))
                        if k == 0:
                            sx = d
                        elif k == 1:
                            sy = d
                        else:
                            sz = d
                if None in (sx, sy, sz):
                    # A 45 deg box has no axis to line up with, so report its
                    # own-frame extents and let rotY45 flag the orientation.
                    sx, sy, sz = dims
                return "BOX", {"pos": transform(c.x, c.y, c.z),
                               "dx": sx, "dy": sy, "dz": sz}

    # --- TUBE: inner and outer cylinder, two flat ends -------------------
    if surf == {"Cylinder", "Plane"} and len(cyls) == 2 and len(solid.Faces) == 4:
        radii = sorted({round(f.Surface.Radius, 6) for f in cyls})
        if len(radii) == 2:
            ax = cyls[0].Surface.Axis
            pts = [v.Point for v in solid.Vertexes]
            pr = [p.x * ax.x + p.y * ax.y + p.z * ax.z for p in pts]
            length = max(pr) - min(pr)
            # Confirm it is a plain hollow cylinder and not, say, a segment:
            # the volume must equal the full annulus. Both tubes here match
            # to 1.00000.
            annulus = math.pi * (radii[1] ** 2 - radii[0] ** 2) * length
            if annulus > 0 and abs(solid.Volume / annulus - 1.0) < TOL_VOL:
                b = solid.BoundBox
                return "TUBE", {"pos": transform(b.Center.x, b.Center.y, b.Center.z),
                                "rmin": radii[0], "rmax": radii[1], "dz": length,
                                # Report the enclosing extents too, so the shape
                                # key and the colour lookup see real sizes.
                                "dx": 2.0 * radii[1], "dy": 2.0 * radii[1],
                                "axis": transform_dir(ax.x, ax.y, ax.z)}

    # --- BOX_HOLE: a block with bores drilled through it -----------------
    # Each cylinder becomes a G4Tubs to subtract. The outer size falls back to
    # the bounding box, which is safe here because these solids are all
    # axis-aligned once de-tilted; a 45 deg one is flagged by rotY45.
    if cyls and "Plane" in surf:
        b = solid.BoundBox
        bores = []
        for f in cyls:
            c = f.Surface.Center
            bores.append({"r": f.Surface.Radius,
                          "pos": transform(c.x, c.y, c.z),
                          "axis": transform_dir(f.Surface.Axis.x,
                                                f.Surface.Axis.y,
                                                f.Surface.Axis.z)})
        return "BOX_HOLE", {"pos": transform(b.Center.x, b.Center.y, b.Center.z),
                            "dx": b.XLength, "dy": b.YLength, "dz": b.ZLength,
                            "bores": bores}

    # --- PRISM / OTHER ---------------------------------------------------
    # Planar but not box-like: wedges, chamfered and skewed blocks. These need
    # G4Trap, G4GenericTrap or G4ExtrudedSolid and are left for a human to
    # decide, so only their envelope is reported. Nothing in this model falls
    # through to OTHER.
    b = solid.BoundBox
    kind = "PRISM" if surf == {"Plane"} else "OTHER"
    return kind, {"pos": transform(b.Center.x, b.Center.y, b.Center.z),
                  "dx": b.XLength, "dy": b.YLength, "dz": b.ZLength}


def shape_key(solid, kind, p):
    """Identity of a block, independent of where it sits or how it is turned.

    Two solids share a key when they are the same part, so the key is what
    makes deduplication work: 224 solids collapse to 45 shapes, and the 146
    copies of the 2x4x8 inch lead brick become one definition however they are
    oriented.

    Sorting the dimensions is what makes it rotation-free -- a block turned on
    its side yields the same triple. The other components keep genuinely
    different parts apart: surface signature, bore radii, and volume. That last
    one matters more than it looks. Several blocks share an outside size but
    differ by a drilled hole (1048772 vs 945809 mm^3); without volume in the
    key they would merge and a bore would be lost.
    """
    sig = tuple(sorted(collections.Counter(
        f.Surface.__class__.__name__ for f in solid.Faces).items()))
    dims = tuple(sorted(round(v, 2) for v in
                        (p.get("dx", 0.0), p.get("dy", 0.0), p.get("dz", 0.0))))
    radii = tuple(sorted({round(f.Surface.Radius, 3) for f in cylinders(solid)}))
    return (kind, sig, dims, radii, round(solid.Volume, 1))


# ----------------------------------------------------------------- material

def _step_solid_colours(path):
    """Colour of every solid in a STEP file, in Part.read order.

    The styling chain STYLED_ITEM -> PRESENTATION_STYLE_ASSIGNMENT -> ... ->
    COLOUR_RGB attaches a colour to each MANIFOLD_SOLID_BREP. Tying those
    entities back to kernel solids is the hard part: the entity ids are not in
    solid order (only 173 of 224 line up), and centroids computed from the
    STEP text are contaminated for about a quarter of solids, because the
    entity graph shares references between neighbouring bodies.

    What does hold: bucket both sides by face count, and within a bucket the
    n-th entity is the n-th solid. Checked on the simplified file, where the
    surface-type signature then agrees for all 224 solids, and on the full
    file, where every bucket matches in size and all 390 solids get a colour.

    Approaches that were tried and did not work, so they are not retried here:
    reading colours through FreeCAD (ImportGui cannot load in a console
    process, and a plain Import exposes no colour property); pairing by entity
    order (173 of 224); pairing by centroids read from the STEP text (51 of 224
    are contaminated, by up to 95 mm); and pairing by dimension class (67 of 83
    classes never match, for the same reason).

    Returns (entities_by_face_count, colour_by_entity_id).
    """
    # Entity table: "#123=PLANE(...)" -> {"123": "PLANE(...)"}. Parsing the
    # text directly, because the kernel discards colour on import.
    txt = open(path, "r", errors="replace").read()
    ent = dict(re.findall(r"#(\d+)=([A-Z_0-9]+\([^;]*\));", txt, re.S))

    def etype(k):
        """Entity type of "#k", e.g. "ADVANCED_FACE"."""
        m = re.match(r"([A-Z_0-9]+)\(", ent.get(k, ""))
        return m.group(1) if m else "?"

    def refs(s):
        """The "#nnn" ids an entity mentions."""
        return re.findall(r"#(\d+)", s)

    # Colour definitions, keyed by id. NX writes names like 'granite gray'.
    colours = {}
    for k, v in ent.items():
        if v.startswith("COLOUR_RGB"):
            m = re.match(r"COLOUR_RGB\('([^']*)'", v)
            if m:
                colours[k] = m.group(1)

    def shell_faces(k):
        """Faces a solid owns, following one link only.

        Deliberately not recursive. Walking deeper picks up geometry shared
        with neighbouring bodies, which is exactly what corrupted the earlier
        attempts. One hop to the shell, then its own faces, stays clean: the
        face-count histogram from this matches the kernel's exactly.
        """
        for c in refs(ent[k]):
            if etype(c) in ("CLOSED_SHELL", "ORIENTED_CLOSED_SHELL"):
                return [f for f in refs(ent[c]) if etype(f) == "ADVANCED_FACE"]
        return []

    # Walk each styled item to the colour at the end of its style chain:
    # STYLED_ITEM -> PRESENTATION_STYLE_ASSIGNMENT -> SURFACE_STYLE_USAGE
    #   -> SURFACE_STYLE_FILL_AREA -> FILL_AREA_STYLE_COLOUR -> COLOUR_RGB.
    # The intermediate steps vary between exporters, so rather than hard-code
    # them, search the subtree for the first colour. Styled items that point
    # at curves or points instead of solids are skipped.
    msb_colour = {}
    for k, v in ent.items():
        if not v.startswith("STYLED_ITEM"):
            continue
        r = refs(v)
        if len(r) < 2 or etype(r[-1]) != "MANIFOLD_SOLID_BREP":
            continue
        seen, stack, col = set(), [r[0]], None
        while stack:
            c = stack.pop()
            if c in seen:
                continue
            seen.add(c)
            if c in colours:
                col = colours[c]
                break
            stack.extend(refs(ent.get(c, "")))
        msb_colour[r[-1]] = col

    # Bucket the solid entities by face count, in declaration order. The
    # ordering within a bucket is what the pairing relies on.
    entities = collections.defaultdict(list)
    for k in sorted((k for k in ent if etype(k) == "MANIFOLD_SOLID_BREP"),
                    key=lambda k: int(k)):
        entities[len(shell_faces(k))].append(k)
    return entities, msb_colour


def solid_colours(path, solids):
    """Colour per kernel solid index, via the face-count bucket pairing.

    Bucket the kernel solids the same way the STEP entities were bucketed,
    then pair them off position by position. A bucket whose two sides differ
    in size is skipped rather than guessed at; that never happens for either
    of these files.
    """
    entities, msb_colour = _step_solid_colours(path)
    buckets = collections.defaultdict(list)
    for i, s in enumerate(solids):
        buckets[len(s.Faces)].append(i)
    out = {}
    for fc, ids in buckets.items():
        ents = entities.get(fc, [])
        if len(ents) != len(ids):
            continue                      # bucket disagrees: leave uncoloured
        for k, i in zip(ents, ids):
            col = msb_colour.get(k)
            if col:
                out[i] = col
    return out


def build_colour_map():
    """Colour for every solid of the simplified model.

    Primary route: match each simplified solid to its counterpart in the
    colour-coded full assembly by bounding-box overlap (222 of 224 agree to
    better than 1%; the two stragglers are parts the simplification reshaped,
    and their nearest neighbour is still unambiguous), then read that solid's
    colour. A dimension-keyed lookup is kept as a fallback for anything the
    overlap match cannot place.

    Both files sit in the same CAD frame, so the comparison needs no transform
    between them. This runs before the Mu2e transform for that reason.

    Returns (colour_by_simplified_index, match_score, fallback_lookup).
    """
    simple = Part.read(SIMPLE)
    full = Part.read(FULL)
    full_colour = solid_colours(FULL, full.Solids)

    # Pre-compute the full assembly's boxes once: this is a 224 x 390 sweep.
    boxes = []
    for j, s in enumerate(full.Solids):
        b = s.BoundBox
        boxes.append((j, b.XMin, b.XMax, b.YMin, b.YMax, b.ZMin, b.ZMax))

    def overlap(a0, a1, b0, b1):
        """Length shared by two intervals, 0 if they miss each other."""
        return max(0.0, min(a1, b1) - max(a0, b0))

    per_solid, scores = {}, {}
    for i, s in enumerate(simple.Solids):
        b = s.BoundBox
        own = max(b.XLength * b.YLength * b.ZLength, 1e-9)
        best, best_j = 0.0, None
        for (j, x0, x1, y0, y1, z0, z1) in boxes:
            iv = (overlap(b.XMin, b.XMax, x0, x1)
                  * overlap(b.YMin, b.YMax, y0, y1)
                  * overlap(b.ZMin, b.ZMax, z0, z1))
            if iv <= 0:
                continue
            # Divide by the LARGER of the two volumes, so a small part sitting
            # inside a big one does not score well. A near-1.0 score therefore
            # means "same size and same place", not merely "contained in".
            other = max((x1 - x0) * (y1 - y0) * (z1 - z0), 1e-9)
            sc = iv / max(own, other)
            if sc > best:
                best, best_j = sc, j
        if best_j is not None:
            col = full_colour.get(best_j)
            if col:
                per_solid[i] = col
                scores[i] = best      # kept, and written out, so weak matches show

    # Fallback lookup, keyed on whole-millimetre sorted dimensions. Built from
    # kernel geometry, not from the STEP text, so it does not inherit the
    # contamination described at _step_solid_colours.
    lut = collections.defaultdict(collections.Counter)
    for j, col in full_colour.items():
        b = full.Solids[j].BoundBox
        key = ",".join(str(int(round(x)))
                       for x in sorted((b.XLength, b.YLength, b.ZLength)))
        lut[key][col] += 1
    return per_solid, scores, {k: c.most_common(1)[0][0] for k, c in lut.items()}


def colour_from_lut(dims, lut):
    """Fallback: colour of the part with these dimensions.

    Only reached when the overlap match found nothing, which does not happen
    for the current files. Tries an exact dimension class first, then the
    nearest one within 3 mm -- enough slack for a fillet the simplification
    removed, too little to reach a different part. Returns ("", "unmatched")
    rather than guessing when nothing is close.
    """
    key = ",".join(str(int(round(x))) for x in sorted(dims))
    if key in lut:
        return lut[key], "dims"
    tgt = sorted(dims)
    best, bestd = None, None
    for k, v in lut.items():
        cand = [float(x) for x in k.split(",")]
        if len(cand) != 3:
            continue
        # Worst single-dimension disagreement, so one badly-off edge is enough
        # to reject a candidate.
        d = max(abs(a - b) for a, b in zip(sorted(cand), tgt))
        if bestd is None or d < bestd:
            best, bestd = v, d
    if best is not None and bestd <= 3.0:
        return best, "dims~%.1fmm" % bestd
    return "", "unmatched"


# --------------------------------------------------------------------- name

def suggest_name(kind, p, colour, seq, count, spread):
    """A starting point for the volume name, following the Offline scheme.

    Offline names shielding as {Material}{Side}wall{N}PV -- CopperLwallPV,
    LeadTwall1PV, BPRwall2PV. That only reads well for a panel that sits on one
    side of the house. Most of this model is stacked brick: 146 copies of one
    2x4x8 inch block scattered over more than a metre in every direction, where
    a side label would be actively misleading. So repeated, widely spread
    shapes are named for what they are, and only shapes that stay put get a
    side. Edit these; they are a starting point, not an authority.
    """
    word = MATERIAL_WORD.get(colour, "Block")
    x, y, z = p["pos"]
    dx, dy, dz = p.get("dx", 0.0), p.get("dy", 0.0), p.get("dz", 0.0)

    if kind == "TUBE":
        return "%sTube%dPV" % (word, seq)

    # Repeated and spread out -> stock brick. Both conditions matter: a shape
    # can repeat a few times and still belong to one wall, and a lone block
    # can sit far from the others.
    if count >= 3 and spread > 200.0:
        inches = sorted(round(v / 25.4) for v in (dx, dy, dz))
        name = "%sBrick%dx%dx%d" % (word, inches[0], inches[1], inches[2])
        # A block of the same outside size but with holes drilled through it is
        # a different part, so say so rather than letting the names collide.
        if p.get("bores"):
            name += "Bored%d" % len(p["bores"])
        return name + "PV"

    # Otherwise call it a panel and guess which wall it belongs to: a slab is
    # thin across the wall it forms, and its sign says which side. Rough by
    # nature -- check these before relying on them. The shape id keeps the
    # names unique even when the guess repeats.
    thin = min(range(3), key=lambda i: (dx, dy, dz)[i])
    if thin == 1:
        side = "T" if y > 0 else "B"          # thin in y -> roof or floor
    elif thin == 0:
        side = "L" if x > 0 else "R"          # thin in x -> side wall
    else:
        side = "F" if z < 300 else "Bk"       # thin in z -> front or back
    return "%s%swall%dPV" % (word, side, seq)


# --------------------------------------------------------------------- main

def main():
    """Read, classify, deduplicate, then write the three CSVs."""
    os.makedirs(OUTDIR, exist_ok=True)

    print("reading %s" % os.path.basename(SIMPLE))
    shape = Part.read(SIMPLE)
    solids = shape.Solids
    print("  %d solids, %d faces" % (len(solids), len(shape.Faces)))

    # NX part numbers (AI-129832-A_1) come from a document import; Part.read
    # alone gives geometry with no names. The two agree solid for solid, which
    # the volume check confirms before the labels are trusted. Names are a
    # convenience, so failure here is a warning rather than an error.
    labels = {}
    try:
        import Import
        doc = FreeCAD.newDocument("stm_labels")
        Import.insert(SIMPLE, "stm_labels")
        feats = [o for o in doc.Objects if o.TypeId == "Part::Feature"]
        if len(feats) == len(solids) and all(
                abs(o.Shape.Volume - s.Volume) < 1e-6
                for o, s in zip(feats, solids)):
            labels = {i: o.Label for i, o in enumerate(feats)}
            print("  matched %d NX part names" % len(labels))
        else:
            print("  WARNING: part names did not align; omitting them")
        FreeCAD.closeDocument(doc.Name)
    except Exception as exc:
        print("  WARNING: could not read part names (%s)" % exc)

    print("recovering materials from the colour-coded full assembly")
    per_solid, scores, lut = build_colour_map()
    print("  %d of %d solids matched to a coloured part" % (len(per_solid), len(solids)))
    print("  %d dimension classes kept as fallback" % len(lut))

    # Classify every solid and file it under its shape key. Solids sharing a
    # key are the same block in different places, so the first one seen
    # supplies the definition and the rest become placements.
    groups = collections.OrderedDict()
    placements = []
    bores = []
    kinds = collections.Counter()

    for i, s in enumerate(solids):
        kind, p = classify(s)
        kinds[kind] += 1
        key = shape_key(s, kind, p)
        if key not in groups:
            groups[key] = {"kind": kind, "p": p, "members": [],
                           "oblique": oblique(s)}
        groups[key]["members"].append(i)
        placements.append((i, key, p, s))

    print("  %d distinct shapes" % len(groups))

    # One row per distinct shape: the block definitions a constructSTM.cc
    # would declare once and place repeatedly.
    order = {k: n for n, k in enumerate(groups)}
    shape_rows = []
    for key, g in groups.items():
        n = order[key]
        p = g["p"]
        dims = (p.get("dx", 0.0), p.get("dy", 0.0), p.get("dz", 0.0))

        # Colour of the group: the members agree, so take the commonest and
        # fall back to the dimension lookup only if none of them matched.
        votes = collections.Counter(per_solid[i] for i in g["members"]
                                    if i in per_solid)
        if votes:
            colour = votes.most_common(1)[0][0]
            how = "part" if len(votes) == 1 else "part(mixed)"
        else:
            colour, how = colour_from_lut(dims, lut)

        # How far apart the copies sit: a brick stacked all over the house
        # should not be named after one wall.
        pts = [placements[i][2]["pos"] for i in g["members"]]
        spread = max(max(q[a] for q in pts) - min(q[a] for q in pts)
                     for a in range(3)) if len(pts) > 1 else 0.0

        name = suggest_name(g["kind"], p, colour, n, len(g["members"]), spread)
        shape_rows.append({
            "shape_id": n,
            "type": g["kind"],
            "count": len(g["members"]),
            "name": name,
            "material": colour,
            "material_hint": COLOUR_HINT.get(colour, ""),
            "material_src": how,
            "dx": round(dims[0], 3),
            "dy": round(dims[1], 3),
            "dz": round(dims[2], 3),
            "rmin": round(p["rmin"], 3) if "rmin" in p else "",
            "rmax": round(p["rmax"], 3) if "rmax" in p else "",
            "nbores": len(p.get("bores", [])) or "",
            "rotY45": "yes" if g["oblique"] else "",
            "in_x": round(dims[0] / 25.4, 3),
            "in_y": round(dims[1] / 25.4, 3),
            "in_z": round(dims[2] / 25.4, 3),
        })

    # One row per solid: where each block goes, in Mu2e coordinates relative
    # to _STMShieldingRef, plus its bores if it has any.
    place_rows = []
    for i, key, p, s in placements:
        x, y, z = p["pos"]
        place_rows.append({
            "solid_id": i,
            "shape_id": order[key],
            "part": labels.get(i, ""),
            "x": round(x, 3),
            "y": round(y, 3),
            "z": round(z, 3),
            "rotY45": "yes" if oblique(s) else "",
            "volume": round(s.Volume, 1),
            # Overlap with the full-assembly part the material came from.
            # Well below 1.0 means simplification reshaped this solid, so the
            # material is a judgement call rather than a direct read.
            "material_match": round(scores[i], 3) if i in scores else "",
        })
        for b in p.get("bores", []):
            bores.append({
                "solid_id": i,
                "shape_id": order[key],
                "r": round(b["r"], 3),
                "x": round(b["pos"][0], 3),
                "y": round(b["pos"][1], 3),
                "z": round(b["pos"][2], 3),
                "ax": round(b["axis"][0], 4),
                "ay": round(b["axis"][1], 4),
                "az": round(b["axis"][2], 4),
            })

    def write(path, rows):
        """Write rows to CSV, taking the column order from the first row."""
        with open(path, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
            w.writeheader()
            w.writerows(rows)
        print("wrote %s (%d rows)" % (path, len(rows)))

    print()
    write(os.path.join(OUTDIR, "stm_shapes.csv"), shape_rows)
    write(os.path.join(OUTDIR, "stm_placements.csv"), place_rows)
    if bores:
        write(os.path.join(OUTDIR, "stm_bores.csv"), bores)

    print()
    print("classification:")
    for k, v in kinds.most_common():
        print("  %-9s %4d" % (k, v))
    print("  %-9s %4d" % ("TOTAL", sum(kinds.values())))

    reused = [r for r in shape_rows if r["count"] > 1]
    print()
    print("shapes defined once and placed many times: %d definitions -> %d solids"
          % (len(reused), sum(r["count"] for r in reused)))
    for r in sorted(reused, key=lambda r: -r["count"])[:6]:
        print("  x%-4d %-16s %8.1f x %6.1f x %6.1f mm  (%s)"
              % (r["count"], r["material"] or "?", r["dx"], r["dy"], r["dz"],
                 r["name"]))

    unmatched = [r for r in shape_rows if r["material_src"] == "unmatched"]
    if unmatched:
        print()
        print("no colour found for %d shapes (material left blank):" % len(unmatched))
        for r in unmatched:
            print("  shape %-3d %8.1f x %6.1f x %6.1f mm  x%d"
                  % (r["shape_id"], r["dx"], r["dy"], r["dz"], r["count"]))


if __name__ == "__main__":
    main()
