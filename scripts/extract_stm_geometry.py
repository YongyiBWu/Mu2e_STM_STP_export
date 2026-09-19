"""
Extract the STM shield-house geometry from the NX STEP export into a
deduplicated, classified inventory in Mu2e coordinates.

Run with FreeCAD's bundled interpreter (FreeCAD must be imported before Part):

    "C:/Program Files/FreeCAD 1.1/bin/python.exe" scripts/extract_stm_geometry.py

Inputs
    STM_STP_files/F10269585_G4_Shield_House_3.stp   (224 solids)
    STM_STP_files/F10258491--_1-Shield House Square.stp          (390 solids)

Outputs (output/)
    stm_shapes.csv      one row per distinct shape: the block definitions
    stm_placements.csv  one row per solid: which shape goes where
    stm_bores.csv       cylindrical cuts, for the G4SubtractionSolid cases
    stm_prisms.csv      cap outlines for the G4ExtrudedSolid cases

Bore positions
    stm_bores.csv gives each bore as dx/dy/dz RELATIVE to its block's centre.
    That is what a G4SubtractionSolid needs, and unlike an absolute position it
    stays valid for every placement of a shape -- the 146 identical bricks share
    one definition, so a world position measured from one of them means nothing
    for the other 145. Absolute x/y/z is carried alongside for cross-checking
    against the CAD only.

    The centre of a bore is NOT a cylindrical face's Surface.Center: OCC puts
    that at the surface's parametric origin, which can sit well outside the
    block. It is recovered by projecting the face's vertices onto its own axis
    and taking the midpoint of the span.

Coordinate transform (CAD -> Mu2e)
    1. De-tilt by 0.041591 deg about z. The NX export carries a spurious
       rotation: planar normals cluster at 0.0416 and 89.9584 deg mod 90, one
       global value, not per-part scatter. It is a CAD error, so it is removed.
       Afterwards 251 of 1364 planar faces are off-axis, and those are the
       genuine 45 deg features.

       Mind the sign: the rotation applied is +TILT_DEG, not -TILT_DEG, because
       step 2 negates x afterwards and that mirror flips the sign of an angle
       in the xy plane. See the comment on _T -- getting this backwards leaves
       2x the tilt in every direction and is nearly invisible downstream.
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
                      "F10269585_G4_Shield_House_3.stp")
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

# Rotating by +TILT_DEG undoes the tilt the exporter baked in -- note the sign.
#
# The obvious choice is -TILT_DEG, and it is wrong, because the rotation is not
# the last thing that happens to a direction: x is negated afterwards (the 180
# deg flip about y). Negating x mirrors the vector, which flips the SIGN of its
# angle in the xy plane. So a direction sitting at +tilt arrives at -tilt from
# the mirror alone, and pre-rotating by -tilt takes it to -2*tilt instead of 0.
#
# The error is invisible to is_axis_aligned(), since 2*tilt displaces a unit
# vector by only 1.05e-6 against a TOL_AXIS of 1e-4, and the BOX paths use
# bounding-box dimensions rather than the vector itself. It shows up only where
# a transformed direction reaches a CSV: bore axes and prism basis vectors.
#
# Checked against all 1364 planar normals: with +TILT_DEG, 251 are left
# off-axis (the genuine 45 deg faces); with -TILT_DEG, 856 are.
_T = math.radians(TILT_DEG)


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


# Caps count as congruent within this relative area difference. Not exact
# equality: solid 178's two caps differ by 1 part in 5e5 -- real CAD noise on a
# face of 4859 mm^2 -- and demanding exactness would reject a genuine
# extrusion and force it back to a bounding envelope.
TOL_CAP_AREA = 1e-4


# Coplanarity tolerance for grouping the faces of one cap, mm. Solid 220's
# cap pieces sit 0.022 mm apart because the solid straddles the piecewise tilt
# -- part modelled tilted, part square -- so exact coplanarity is too strict.
# Still far below the smallest real feature (12.7 mm).
TOL_PLANE = 0.05

# Volume tolerance for a SPLIT cap only; single-face caps keep TOL_VOL. The
# step that splits a cap also stops the solid being an exact sweep: solid 220
# reads 1.4968e-4 against area*length. See the note at the check itself.
TOL_VOL_SPLIT = 5e-4


def _outline_area(face):
    """Area of a face's outer boundary, ignoring holes drilled through it.

    face.Area is the material left AFTER a bore breaks the surface, so two
    congruent caps report different areas when a hole pierces one of them.
    Adding the inner wires back compares the outlines themselves.
    """
    inner = 0.0
    for w in face.Wires:
        if not w.isSame(face.OuterWire):
            try:
                inner += Part.Face(w).Area
            except Exception:
                pass
    return face.Area + inner


def _union_boundary(faces, nd=3):
    """Ordered boundary points of several coplanar faces, or None.

    Edges interior to the union appear twice (once from each face) and edges on
    the boundary appear once, so dropping the doubled ones and walking what is
    left gives the outline in true order. Ordering by angle about the centroid
    would be simpler and is wrong: 5 of the 13 caps here are non-convex.

    Returns the ORIGINAL FreeCAD points, full precision. Rounding is used only
    to decide which vertices are the same one: returning the rounded values
    instead silently truncated stored outlines (526.5575 became 526.5571).
    """
    def key(p):
        return (round(p.x, nd), round(p.y, nd), round(p.z, nd))

    exact = {}
    seen = collections.defaultdict(list)
    for f in faces:
        for e in f.OuterWire.Edges:
            a, b = e.Vertexes[0].Point, e.Vertexes[-1].Point
            ka, kb = key(a), key(b)
            exact.setdefault(ka, a)
            exact.setdefault(kb, b)
            seen[frozenset((ka, kb))].append((ka, kb))
    bnd = [v[0] for v in seen.values() if len(v) == 1]
    if not bnd:
        return None

    adj = collections.defaultdict(list)
    for a, b in bnd:
        adj[a].append(b)
        adj[b].append(a)
    # A clean ring: every corner joins exactly two boundary edges.
    if any(len(v) != 2 for v in adj.values()):
        return None

    start = bnd[0][0]
    order = [start]
    prev, cur = None, start
    while True:
        nxt = [n for n in adj[cur] if n != prev]
        if not nxt:
            return None
        step = nxt[0]
        if step == start:
            break
        order.append(step)
        prev, cur = cur, step
        if len(order) > len(adj):
            return None
    return [exact[k] for k in order]


def bore_volume(solid):
    """Material removed by this solid's cylindrical faces, mm^3.

    Each bore's length is its own face's vertex span along its axis, the same
    measurement the BOX_HOLE branch uses for depth, so a blind hole is not
    charged the full thickness of the block.
    """
    total = 0.0
    for f in cylinders(solid):
        ax = f.Surface.Axis
        pts = [v.Point for v in f.Vertexes]
        if not pts:
            continue
        pr = [p.x * ax.x + p.y * ax.y + p.z * ax.z for p in pts]
        total += math.pi * f.Surface.Radius ** 2 * (max(pr) - min(pr))
    return total


def extrusion(solid, planar_only=False, ref_volume=None):
    """Cap polygon and length if this solid is an extrusion, else None.

    planar_only=True ignores cylindrical faces, for a block that is an
    extrusion once its bores are set aside. ref_volume then supplies the
    volume the sweep is checked against -- solid.Volume plus what the bores
    removed -- because a drilled solid is lighter than its own outline sweeps.

    An extrusion has one axis along which exactly two faces are perpendicular
    (the caps, congruent and parallel) and every remaining face is parallel
    (the walls). That is precisely what G4ExtrudedSolid takes: a 2D outline
    swept a given distance.

    Returns the cap outline as 2D points in the cap's own plane, with the two
    in-plane basis vectors and the axis, all in Mu2e coordinates, so a
    constructSTM.cc can build the solid and orient it without going back to
    the STEP file.

    The volume is checked against area*length: a solid can pass the face test
    and still be a sheared prism whose volume falls short.
    """
    faces = []
    for f in solid.Faces:
        if f.Surface.__class__.__name__ != "Plane":
            if planar_only:
                continue
            return None
        n = f.normalAt(0, 0)
        n.normalize()
        faces.append((f, n))
    vol = solid.Volume if ref_volume is None else ref_volume

    for _, cand in faces:
        caps = [f for f, n in faces if abs(abs(n.dot(cand)) - 1.0) < TOL_UNIT]
        walls = [f for f, n in faces if abs(n.dot(cand)) < TOL_UNIT]
        if len(caps) < 2 or len(caps) + len(walls) != len(faces):
            continue

        # One cap can be SPLIT across several faces, so group the cap-normal
        # faces by where they sit along the candidate axis instead of demanding
        # exactly two of them.
        #
        # Solid 220 is the case: swept along y, both of its caps are two faces,
        # which the old "exactly two" test read as four caps and rejected. Its
        # two pieces are also 0.022 mm out of plane, because part of that solid
        # was modelled tilted and part square -- it straddles the piecewise
        # tilt -- so TOL_PLANE has to absorb a real step, not just noise. It
        # stays far below the smallest true feature here (12.7 mm).
        groups = []
        for f in caps:
            p = f.Vertexes[0].Point
            off = p.x * cand.x + p.y * cand.y + p.z * cand.z
            for g in groups:
                if abs(g[0] - off) < TOL_PLANE:
                    g[1].append(f)
                    break
            else:
                groups.append((off, [f]))
        if len(groups) != 2:
            continue

        # Congruence on the OUTLINE, adding back what a bore removed: a cap a
        # hole breaks through reports less area than its twin, and solid 220's
        # z-faces differ by exactly one bore's 4053.7 mm^2 that way.
        near, far = groups[0][1], groups[1][1]
        a0 = sum(_outline_area(f) for f in near)
        a1 = sum(_outline_area(f) for f in far)
        if abs(a0 - a1) / max(a0, a1) > TOL_CAP_AREA:
            continue

        # The outline, as points to walk.
        #
        # A single-face cap keeps its own OuterWire verbatim -- same order,
        # same start vertex, same full precision. That is deliberate: routing
        # every prism through the stitching path instead changed 12 of 13
        # outlines (flipped bases, moved anchors, truncated coordinates) while
        # still not recovering 220. Stitching is only for the case that needs
        # it, so everything that already worked is untouched by construction.
        if len(near) == 1:
            cap_pts = [v.Point for v in near[0].OuterWire.OrderedVertexes]
        else:
            cap_pts = _union_boundary(near)
            if cap_pts is None:
                continue

        pts = [v.Point for v in solid.Vertexes]
        pr = [p.x * cand.x + p.y * cand.y + p.z * cand.z for p in pts]
        length = max(pr) - min(pr)

        # A split cap gets a looser volume bound, because the step that splits
        # it makes the solid not quite a sweep.
        #
        # Solid 220's cap pieces sit 0.022 mm out of plane -- it straddles the
        # piecewise tilt -- and over its 3932.5 mm^2 profile that wedge is 81.06
        # mm^3, so area*length overstates the true volume by 1.4968e-4. That is
        # real geometry, not a bad measurement: bore_volume() was checked
        # against the material actually removed (holes plugged and re-measured)
        # and agreed to -0.000 mm^3.
        #
        # Deliberately NOT applied when both caps are single faces. TOL_VOL is
        # what stops a sheared prism from passing as a clean one, and every
        # solid that already classifies keeps it untouched; only the case whose
        # cause is identified is allowed the slack, and 5e-4 still rejects any
        # shear worth the name.
        tol_vol = TOL_VOL_SPLIT if (len(near) > 1 or len(far) > 1) else TOL_VOL
        if length <= 0 or abs(a0 * length / vol - 1.0) > tol_vol:
            continue

        # Sweep midpoint along the axis: the plane the origin sits on.
        mid = (max(pr) + min(pr)) / 2.0

        # In-plane axes for expressing the cap outline in 2D.
        #
        # Any pair square to the extrusion axis is geometrically valid, but the
        # CHOICE decides how easy the later rotation is to reason about. Prefer
        # an edge that already runs along a world axis once transformed: then
        # local x (or y) coincides with Mu2e x (or y) and the rotation is a
        # plain permutation instead of an arbitrary in-plane spin.
        #
        # Deliberately NOT the longest edge -- that rule fails here. Shape 23's
        # longest edge is 744.665 mm and already on x, so it would change
        # nothing, while 6 of 11 prisms have their longest edge at 45 or 90 deg
        # to everything.
        # Basis candidates come from the FACE's own edges when the cap is a
        # single face, exactly as before. Walking cap_pts instead reverses the
        # traversal direction on some caps, which flips the sign of u (and so
        # of w and the whole outline) -- it silently changed 6 of 13 prisms
        # that way. Only a stitched cap, which has no single face to ask, uses
        # the walked boundary.
        edges = []
        if len(near) == 1:
            for e in near[0].Edges:
                d = e.Vertexes[-1].Point.sub(e.Vertexes[0].Point)
                if d.Length > 1e-9:
                    d.normalize()
                    if abs(d.dot(cand)) < TOL_UNIT:
                        edges.append(d)
        else:
            for i in range(len(cap_pts)):
                d = cap_pts[(i + 1) % len(cap_pts)].sub(cap_pts[i])
                if d.Length > 1e-9:
                    d.normalize()
                    if abs(d.dot(cand)) < TOL_UNIT:
                        edges.append(d)
        if not edges:
            return None

        # Rank: an edge whose TRANSFORMED direction lies on a world axis wins;
        # x is preferred over y so the choice is deterministic when both exist.
        def rank(d):
            t = transform_dir(d.x, d.y, d.z)
            m = max(range(3), key=lambda i: abs(t[i]))
            if abs(abs(t[m]) - 1.0) > TOL_AXIS:
                return (2, 0)              # not axis-aligned at all
            return (0 if m == 0 else 1, m)  # x best, then y/z

        u = min(edges, key=rank)
        w = cand.cross(u)
        w.normalize()
        # If w landed on a world axis but u did not, swap them so the aligned
        # edge is local x rather than local y.
        if rank(u)[0] == 2 and rank(w)[0] < 2:
            u, w = w, cand.cross(w)
            w.normalize()

        # Walk the cap's outer wire so the polygon comes out ordered, which is
        # what G4ExtrudedSolid needs -- an unordered point set would build a
        # self-intersecting face.
        #
        # The outline is expressed RELATIVE TO THE SOLID'S CENTRE, not as raw
        # projected CAD coordinates. transform() is rotate -> negate x,z ->
        # subtract the anchor; the rotation and negation are linear and so
        # survive being folded into the stored u/w basis, but the anchor
        # translation does not. Projecting untransformed points onto the basis
        # silently dropped it, leaving every outline offset by a constant ~512mm
        # and unplaceable. Subtracting the centre first removes the translation
        # from both sides, so the outline composes with a placement's x,y,z the
        # same way bore offsets do.
        raw = []
        for p in cap_pts:
            raw.append((p.x * u.x + p.y * u.y + p.z * u.z,
                        p.x * w.x + p.y * w.y + p.z * w.z))

        # Origin on a VERTEX, not on the centroid.
        #
        # The outline still has to be expressed relative to something on the
        # solid rather than as raw projected CAD coordinates: transform() is
        # rotate -> negate x,z -> subtract the anchor, and while the rotation
        # and negation fold into the stored u/w basis, the anchor translation
        # does not. Projecting untransformed points left every outline offset
        # by a constant ~512 mm and unplaceable.
        #
        # The centre satisfied that but put (0,0) inside the polygon, which is
        # awkward to reason about when composing rotations by hand. Anchoring
        # on the first ordered vertex keeps the translation removed AND puts
        # the origin on a corner you can point at.
        ou, ow = raw[0]
        outline = [(round(a - ou, 4), round(b - ow, 4)) for a, b in raw]

        # The placement has to move WITH the origin. pos is the point the
        # outline is expressed relative to, so anchoring on a vertex while pos
        # still pointed at the solid's centre put the vertex where the centre
        # belonged -- every prism landed 350mm out and the round-trip overlap
        # fell to 0.0008. The anchor vertex in CAD space is the first ordered
        # vertex projected back onto the basis, plus its component along the
        # sweep axis at the mid-plane, so the origin sits on the cap's vertex
        # at the sweep half-length.
        anchor_cad = (u * ou) + (w * ow) + (cand * mid)
        anchor = transform(anchor_cad.x, anchor_cad.y, anchor_cad.z)
        return {
            "prism_axis": transform_dir(cand.x, cand.y, cand.z),
            "prism_u": transform_dir(u.x, u.y, u.z),
            "prism_w": transform_dir(w.x, w.y, w.z),
            "prism_len": length,
            "prism_outline": outline,
            "prism_area": a0,
            # Where the outline's (0,0) actually sits, in Mu2e. The caller
            # overrides pos with this so the placement and the origin describe
            # one point.
            "prism_anchor": anchor,
        }
    return None


def canonical_frame(axes, dims, longest_x, symmetric=False):
    """Put a solid in a canonical frame, and give the rotation back to world.

    The shape row and the rotation column have to share one convention or they
    cannot be composed. Deriving them separately is what broke earlier: the
    dimensions were labelled by which WORLD axis each own-axis pointed along,
    while the rotation was built from the solid's own axis ORDER. Those are
    different orderings, so applying the rotation to the dimensions
    double-counted the orientation.

    Here the canonical frame is defined first and the rotation is defined AS
    the map from it to the world, so the two agree by construction.

    The convention, per the Geant4 build being the destination:

        longest_x=False (a plain box)
            the two longest extents lie on x and y, the shortest on z, so a
            slab sits flat and its height is z.
        longest_x=True (anything with a bore, a cylinder, or a swept profile)
            the longest extent lies on x, since those parts are built along
            their axis.

    Right-handedness: a permutation of three axes can be improper (det = -1),
    which would mirror the solid rather than rotate it. If the chosen ordering
    comes out improper, one axis is negated -- that is a rotation of the same
    frame, and the extents are unsigned, so nothing else changes.

    Returns (dims_canonical, rotation) where rotation[row][col] maps the
    canonical frame to Mu2e axes, or (dims, None) when the solid has no
    axis-aligned orientation to express.
    """
    order = sorted(range(3), key=lambda i: -dims[i])       # longest first
    if longest_x:
        # longest -> x, then the remaining two longest -> y, z.
        pick = [order[0], order[1], order[2]]
    else:
        # two longest -> x, y; shortest -> z.
        pick = [order[0], order[1], order[2]]
    cdims = [dims[i] for i in pick]
    cax = [axes[i] for i in pick]

    # Where does each canonical axis point in Mu2e?
    #
    # Two cases. When every canonical axis lands on a world axis the rotation
    # is a signed permutation, and the sign canonicalisation below can use the
    # box's own symmetry. When one does not -- a 45 deg part -- the rotation is
    # still a perfectly ordinary proper rotation, just not a permutation, so
    # emit the exact direction cosines instead of refusing.
    #
    # Refusing was the old behaviour, and it was over-cautious: it made sense
    # while the only consumer was the permutation-based orientation string, but
    # r11..r33 are floats and hold any rotation. Returning None there cost six
    # solids their placement -- three 45 deg lead bricks (171/172/176) and
    # three bored walls -- which then could not be rebuilt from the CSVs, were
    # skipped by the verifier, and reached the G4 emitter as "set by hand".
    dirs = [transform_dir(u.x, u.y, u.z) for u in cax]
    world = []
    for n in dirs:
        if not is_axis_aligned(n):
            world = None
            break
        k = max(range(3), key=lambda i: abs(n[i]))
        world.append((k, 1 if n[k] > 0 else -1))
    if world is not None and len({k for k, _ in world}) != 3:
        world = None

    if world is None:
        # General rotation: column `own` is where that canonical axis points.
        # Orthonormal because the axes came from orthogonal_frame(), which
        # already checked they are mutually square.
        r = [[dirs[own][row] for own in range(3)] for row in range(3)]
        det = (r[0][0] * (r[1][1] * r[2][2] - r[1][2] * r[2][1])
               - r[0][1] * (r[1][0] * r[2][2] - r[1][2] * r[2][0])
               + r[0][2] * (r[1][0] * r[2][1] - r[1][1] * r[2][0]))
        if det < 0:
            # Improper: the measured axes form a left-handed set. Negating one
            # column re-orients the frame without moving any face, exactly as
            # in the permutation case. Extents are unsigned, so cdims stands.
            for row in range(3):
                r[row][2] = -r[row][2]
        # No sign canonicalisation here: it relies on a 180 deg flip about an
        # own axis being a self-map, which is a statement about the signed
        # permutation, not about a general rotation.
        return (cdims, r)

    r = [[0, 0, 0], [0, 0, 0], [0, 0, 0]]
    for own, (k, sign) in enumerate(world):
        r[k][own] = sign

    # det must be +1; if the permutation is improper, flip the z column, which
    # re-orients the canonical frame without changing any measured extent.
    det = (r[0][0] * (r[1][1] * r[2][2] - r[1][2] * r[2][1])
           - r[0][1] * (r[1][0] * r[2][2] - r[1][2] * r[2][0])
           + r[0][2] * (r[1][0] * r[2][1] - r[1][1] * r[2][0]))
    if det < 0:
        for row in range(3):
            r[row][2] = -r[row][2]

    # Canonicalise the signs, for a solid that is symmetric under a 180 deg
    # flip about each of its own axes -- a PLAIN BOX and nothing else here.
    #
    # A rectangular box with three distinct edges has 24 rotational symmetries
    # but only 6 distinct ways of assigning its edges to the world axes: the 6
    # permutations. The other 4 per permutation are 180 deg flips that map the
    # box onto itself, placing identical material and differing only in which
    # end of each edge points which way. Recording the raw signed permutation
    # over-reports -- shape 1 showed 13 distinct orientation strings for 5
    # distinct placements plus 3 unrotatable copies.
    #
    # Flips must come in PAIRS: negating a single axis gives det = -1, a
    # reflection, which would mirror the solid rather than rotate it. Flipping
    # two at a time preserves det and keeps the frame proper, so the normal
    # form is "as few negatives as parity allows" -- 0 or 1.
    #
    # NOT applied to BOX_HOLE: the bores break the symmetry. Flipping shape 9
    # about y or z moves its holes, and a self-overlap test puts that at 0.889,
    # not 1.0. Prisms do not use this path and are not symmetric either
    # (shape 19 self-overlaps at 0.000 about y).
    if symmetric:
        # Reduce to a DETERMINISTIC representative, not merely to "few
        # negatives". Flipping pairs until fewer than two remain leaves 0 or 1,
        # but which axis keeps the surviving negative depends on the order the
        # pairs were consumed -- so two copies sharing one permutation could
        # normalise differently and still look like distinct orientations.
        # Shape 1 stalled at 8 strings that way, with "x->-x y->z z->y" and
        # "x->x y->z z->-y" both surviving as the same permutation.
        #
        # Parity is the only invariant a pair-flip preserves, so it fixes how
        # many negatives must remain: an even count reduces to 0, an odd count
        # to exactly 1. Put that lone negative on the lowest-indexed axis every
        # time and the representative is unique.
        def flip(own):
            r[world[own][0]][own] = -r[world[own][0]][own]

        neg = [own for own in range(3) if r[world[own][0]][own] < 0]

        # Parity is the only invariant a pair-flip preserves, so it fixes how
        # many negatives survive: an even count reduces to 0, an odd count to
        # exactly 1. The subtlety is WHICH axis keeps that lone negative.
        #
        # Keying on the own-axis index does not work: copies sharing one
        # permutation map their own axes to different world axes, so "lowest
        # own axis" picks differently for each and the representative is not
        # unique. Shape 1 split 55/2 that way, one group holding the negative
        # on own-x and the other on own-z.
        #
        # The world axis is what the copies of a permutation agree on, so
        # choose by that: move the lone negative onto whichever own-axis maps
        # to the LOWEST world axis. Moving it is itself a pair flip (clear it
        # here, set it there), which keeps det = +1.
        if len(neg) % 2 == 0:
            for own in neg:
                flip(own)
        else:
            target = min(range(3), key=lambda own: world[own][0])
            for own in neg:
                if own != target:
                    flip(own)
            if r[world[target][0]][target] > 0:
                flip(target)
    return (cdims, r)


def orientation_label(rot):
    """Readable form of a rotation, e.g. "x->y y->-z z->x".

    Picks the DOMINANT entry in each column, not the first non-zero one.
    A box rotation is a signed permutation, so any non-zero entry is the right
    one and either rule works. A prism's rotation is not: its columns are the
    cap basis and sweep axis, which carry the de-tilt residue as off-axis
    components around 7e-4. "First non-zero" then latches onto that noise and
    reports impossible labels -- every prism came out as "x->x y->-x z->-x",
    three axes all mapped to x, while the matrices themselves were exact
    (det = +1, orthonormal to 1e-13).
    """
    if not rot:
        return ""
    names = "xyz"
    out = []
    for own in range(3):
        col = [rot[world][own] for world in range(3)]
        order = sorted(range(3), key=lambda i: -abs(col[i]))
        best, second = order[0], order[1]
        # A 45 deg part has no permutation label: its axis lies BETWEEN two
        # world axes, both components ~0.707, and picking the larger is then
        # arbitrary -- it produced "x->x y->-x z->y" for solid 135, mapping two
        # canonical axes onto the same world axis, which no rotation can do.
        # Say so instead of inventing a permutation. r11..r33 still carry the
        # exact rotation (det +1, orthonormal); only this readable form is
        # undefined.
        if abs(col[second]) > 0.5 * abs(col[best]):
            return "oblique-45"
        out.append("%s->%s%s" % (names[own],
                                 "-" if col[best] < 0 else "",
                                 names[best]))
    return " ".join(out)


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
        PRISM     all-planar but not a box: wedges, chamfered blocks.
                  Every one in this model is an extrusion, so the cap outline
                  and sweep length are recovered -> G4ExtrudedSolid. See
                  extrusion() and stm_prisms.csv.
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
                # Canonical frame: a plain box lies on its two longest sides
                # with the shortest as its height in z. The rotation back to
                # world comes from the same call, so dimensions and rotation
                # cannot disagree.
                # symmetric=True: a plain box maps onto itself under a 180 deg
                # flip about any of its own axes, so the sign of each axis is
                # not observable and is normalised away. BOX_HOLE below passes
                # False -- its bores make the flips real.
                cdims, rot = canonical_frame(axes, dims, longest_x=False,
                                             symmetric=True)
                return "BOX", {"pos": transform(c.x, c.y, c.z),
                               "dx": cdims[0], "dy": cdims[1], "dz": cdims[2],
                               "rot": rot}

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
                axis = transform_dir(ax.x, ax.y, ax.z)
                # A G4Tubs is built along ITS OWN z, so the placement rotation
                # is whatever carries z onto the measured axis. Both tubes here
                # lie exactly on world z (0.0000 deg), so that rotation is the
                # identity -- but it was never recorded, and a blank r11..r33
                # reached the emitter as "TODO: set by hand" on a solid whose
                # orientation was never in doubt.
                #
                # Derived from the axis rather than assumed, so a tube mounted
                # off-axis gets a real matrix instead of a false identity.
                # Columns are where the tube's own x, y, z land in Mu2e.
                zc = FreeCAD.Vector(*axis)
                if zc.Length < 1e-9:
                    rot = None
                else:
                    zc.normalize()
                    # Any vector not parallel to the axis seeds the transverse
                    # pair; which one is arbitrary, because a tube is a body of
                    # revolution and spinning it about its own axis places the
                    # same material.
                    seed = FreeCAD.Vector(1, 0, 0)
                    if abs(seed.dot(zc)) > 0.9:
                        seed = FreeCAD.Vector(0, 1, 0)
                    xc = seed.sub(zc * seed.dot(zc))
                    xc.normalize()
                    yc = zc.cross(xc)          # right-handed: det = +1
                    rot = [[xc[k], yc[k], zc[k]] for k in range(3)]
                return "TUBE", {"pos": transform(b.Center.x, b.Center.y, b.Center.z),
                                "rmin": radii[0], "rmax": radii[1], "dz": length,
                                # Report the enclosing extents too, so the shape
                                # key and the colour lookup see real sizes.
                                "dx": 2.0 * radii[1], "dy": 2.0 * radii[1],
                                "rot": rot,
                                "axis": axis}

    # --- PRISM_HOLE: a bored block whose outline is not a rectangle ------
    #
    # Caught BEFORE BOX_HOLE, and only when orthogonal_frame() has already
    # failed, so nothing that BOX_HOLE handles correctly is diverted here.
    #
    # These are the solids BOX_HOLE was describing with its bounding box: a
    # block with one end cut at 45 deg, drilled through. Calling that a cuboid
    # overstated solid 162 by 12.66% of its volume, and since no orthogonal
    # frame exists it also got no rotation -- so the CSV could not place it and
    # the verifier skipped it. Swept outline plus bores is what it actually is,
    # and G4ExtrudedSolid minus G4Tubs is how Geant4 builds it.
    if cyls and "Plane" in surf and orthogonal_frame(solid) is None:
        ext = extrusion(solid, planar_only=True,
                        ref_volume=solid.Volume + bore_volume(solid))
        if ext is not None:
            u, w, ax = ext["prism_u"], ext["prism_w"], ext["prism_axis"]
            # Same frame convention as PRISM: cap in x,y and sweep along z,
            # which is the frame G4ExtrudedSolid builds in.
            rot = [[u[k], w[k], ax[k]] for k in range(3)]
            anchor = ext["prism_anchor"]
            us = [p[0] for p in ext["prism_outline"]]
            vs = [p[1] for p in ext["prism_outline"]]
            params = dict(ext)
            params["pos"] = anchor
            params["rot"] = rot
            params["dx"] = max(us) - min(us)
            params["dy"] = max(vs) - min(vs)
            params["dz"] = ext["prism_len"]
            # Bores, in the prism's own frame. The offset is measured from the
            # OUTLINE ANCHOR, not from the solid's centre, because that anchor
            # is what pos names and what the outline's (0,0) sits on -- the
            # same rule the PRISM branch follows. Measuring from the centre
            # here would displace every hole by the anchor-to-centre vector.
            bores = []
            for f in cyls:
                fax = f.Surface.Axis
                pts = [v.Point for v in f.Vertexes]
                if not pts:
                    continue
                pr = [p.x * fax.x + p.y * fax.y + p.z * fax.z for p in pts]
                mid = (max(pr) + min(pr)) / 2.0
                depth = max(pr) - min(pr)
                c = f.Surface.Center
                off = mid - (c.x * fax.x + c.y * fax.y + c.z * fax.z)
                world = transform(c.x + fax.x * off, c.y + fax.y * off,
                                  c.z + fax.z * off)
                rel_w = tuple(world[k] - anchor[k] for k in range(3))
                axis_w = transform_dir(fax.x, fax.y, fax.z)
                # R maps the prism frame to world and is orthonormal, so its
                # transpose takes the world offset back into that frame.
                rel = tuple(sum(rot[k][j] * rel_w[k] for k in range(3))
                            for j in range(3))
                axis_c = tuple(sum(rot[k][j] * axis_w[k] for k in range(3))
                               for j in range(3))
                bores.append({"r": f.Surface.Radius, "pos": world,
                              "rel": rel, "depth": depth, "axis": axis_c})
            params["bores"] = bores
            return "PRISM_HOLE", params

    # --- BOX_HOLE: a block with bores drilled through it -----------------
    # Each cylinder becomes a G4Tubs to subtract. The outer size falls back to
    # the bounding box, which is safe here because these solids are all
    # axis-aligned once de-tilted; a 45 deg one is flagged by rotY45.
    if cyls and "Plane" in surf:
        b = solid.BoundBox
        centre = transform(b.Center.x, b.Center.y, b.Center.z)

        # Fix the canonical frame BEFORE measuring the bores: each bore's
        # offset is expressed in that frame, so the rotation has to exist
        # first. Bored blocks are built along their bore axis, so the longest
        # extent goes on x. A block with no clean orthogonal frame keeps its
        # bounding-box extents and gets no rotation, which is reported rather
        # than guessed at.
        _axes = orthogonal_frame(solid)
        if _axes is not None:
            _dims, _mids = extents_along(solid, _axes)
            _cdims, _rot = canonical_frame(_axes, _dims, longest_x=True)
        else:
            _cdims, _rot = ([b.XLength, b.YLength, b.ZLength], None)

        bores = []
        for f in cyls:
            # A cylindrical face's Surface.Center is a point on the axis, NOT
            # the centre of the drilled hole: OCC puts it at the surface's
            # parametric origin, which can sit far outside the block. Taking it
            # as the bore position put 19 of 32 bores outside their own solid.
            #
            # Project the face's own vertices onto its axis instead and take the
            # midpoint of that span: that is the centre of the actual hole, and
            # it lies on the block by construction.
            ax = f.Surface.Axis
            pts = [v.Point for v in f.Vertexes]
            if pts:
                pr = [p.x * ax.x + p.y * ax.y + p.z * ax.z for p in pts]
                mid = (max(pr) + min(pr)) / 2.0
                depth = max(pr) - min(pr)
            else:
                # A full cylinder can have no vertices; fall back to the
                # bounding box centre projected onto the axis.
                fb = f.BoundBox
                mid = (fb.Center.x * ax.x + fb.Center.y * ax.y
                       + fb.Center.z * ax.z)
                depth = 0.0
            c = f.Surface.Center
            # Slide the surface origin along the axis to the hole's midpoint.
            off = mid - (c.x * ax.x + c.y * ax.y + c.z * ax.z)
            cx = c.x + ax.x * off
            cy = c.y + ax.y * off
            cz = c.z + ax.z * off
            world = transform(cx, cy, cz)
            # Position relative to the block centre, expressed in the SHAPE'S
            # CANONICAL FRAME -- not in world axes.
            #
            # A bore belongs to the shape, not to the world: the 146 identical
            # bricks share one definition, so an absolute position is
            # meaningless for every placement but the one it was measured
            # from. That much was already true. What changed is that dx/dy/dz
            # are now canonical (longest extent on x for a bored block), so a
            # world-frame offset no longer indexes the box it is subtracted
            # from: 11 of 32 bores landed outside their own block.
            #
            # Rotating the offset by R-transpose (R maps canonical -> world,
            # and R is orthonormal so its transpose is its inverse) puts the
            # bore back in the frame the box is built in. The axis goes through
            # the same rotation, or a bore would be drilled along the wrong
            # edge once the block is turned.
            rel_w = (world[0] - centre[0], world[1] - centre[1],
                     world[2] - centre[2])
            axis_w = transform_dir(ax.x, ax.y, ax.z)
            if _rot is not None:
                rel = tuple(sum(_rot[k][j] * rel_w[k] for k in range(3))
                            for j in range(3))
                axis_c = tuple(sum(_rot[k][j] * axis_w[k] for k in range(3))
                               for j in range(3))
            else:
                rel, axis_c = rel_w, axis_w
            bores.append({"r": f.Surface.Radius,
                          "pos": world,
                          "rel": rel,
                          "depth": depth,
                          "axis": axis_c})
        # dx/dy/dz ARE the canonical extents -- not the bounding box. Carrying
        # both invites the shape row and the rotation to disagree, which is the
        # failure this rework exists to remove. The bbox is still available to
        # anything that wants it via the STEP itself.
        return "BOX_HOLE", {"pos": centre, "rot": _rot,
                            "dx": _cdims[0], "dy": _cdims[1], "dz": _cdims[2],
                            "bores": bores}

    # --- PRISM / OTHER ---------------------------------------------------
    # Planar but not box-like: wedges and chamfered blocks. Every one of these
    # in the model is an EXTRUSION -- two congruent parallel caps joined by
    # walls parallel to the axis -- which is exactly G4ExtrudedSolid's model, so
    # the cap polygon and the extrusion length are recovered here rather than
    # leaving a human to re-derive them from the CAD.
    b = solid.BoundBox
    kind = "PRISM" if surf == {"Plane"} else "OTHER"
    params = {"pos": transform(b.Center.x, b.Center.y, b.Center.z),
              "dx": b.XLength, "dy": b.YLength, "dz": b.ZLength}
    if kind == "PRISM":
        ext = extrusion(solid)
        if ext is not None:
            params.update(ext)
            # A prism is NOT oblique just because its faces are not world-
            # aligned. extrusion() already recovers a complete orthonormal
            # right-handed frame -- the cap's in-plane basis u, w and the sweep
            # axis -- verified for all 11 prisms: |u|=|w|=|axis|=1, mutual dots
            # below 1e-13, det = +1 exactly.
            #
            # That frame IS the canonical one, and it matches the convention
            # G4ExtrudedSolid needs: the cap polygon lies in x,y and the sweep
            # runs along z. So u->x, w->y, axis->z, and the rotation taking
            # that frame to Mu2e axes is simply those three vectors as COLUMNS
            # (column j is where canonical axis j lands).
            #
            # Without this every prism was reported oblique and emitted with a
            # nullptr rotation, when its orientation was fully determined all
            # along.
            u, w, ax = ext["prism_u"], ext["prism_w"], ext["prism_axis"]
            params["rot"] = [[u[k], w[k], ax[k]] for k in range(3)]
            # The outline is anchored on a vertex, so the placement must name
            # that same point. Leaving pos at the bounding-box centre put every
            # prism ~350 mm from where it belongs.
            params["pos"] = ext["prism_anchor"]
            # Canonical extents: the cap's span in u and v, and the sweep.
            us = [p[0] for p in ext["prism_outline"]]
            vs = [p[1] for p in ext["prism_outline"]]
            params["dx"] = max(us) - min(us)
            params["dy"] = max(vs) - min(vs)
            params["dz"] = ext["prism_len"]
    return kind, params


# Shape identity carried over from iteration 2, when one is available.
#
# Iteration 3 shrank most parts by tenths of a millimetre -- within
# manufacturing tolerance, but enough that shape_key() no longer merges
# nominally identical blocks: the 2x4x8 lead brick split into 28 "shapes" of
# 203.2 x 101.6 x 50.573 / 50.584 / 50.657 / 50.767 ... and the model went from
# 46 shapes to 86. Solid ids are stable between iterations, so the honest fix
# is to take IDENTITY from v2 and GEOMETRY from v3 rather than invent a looser
# rounding that would merge genuinely different parts elsewhere.
V2_MAP_PATH = os.path.join(OUTDIR, "v2_shape_map.json")


def _load_v2_map():
    try:
        with open(V2_MAP_PATH) as fh:
            return json.load(fh)
    except (IOError, OSError, ValueError):
        return {}


V2_MAP = _load_v2_map()


def _load_solid_materials():
    """solid_id -> Geant4 material, from the material pass.

    Needed for naming: a label that names no material ("Shelf Changed") or
    starts with a digit ("248 Lead Brick") gets the material word prepended,
    and the material is not derivable from the CAD colour -- "medium steel"
    alone covers Al6061, BP, Polyethylene and MildSteel.
    """
    try:
        with open(os.path.join(OUTDIR, "solid_materials.json")) as fh:
            return json.load(fh)
    except (IOError, OSError, ValueError):
        return {}


SOLID_MATERIALS = _load_solid_materials()


def v2_key(solid_id):
    """('v2', shape) for a solid iteration 2 already grouped, else None."""
    row = V2_MAP.get(str(solid_id))
    if not row:
        return None
    return ("v2", row["v2_shape"])


def v2_dims(solid_id):
    """Iteration 2's dimensions for this solid's shape, or None."""
    return V2_MAP.get(str(solid_id))


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

    # WHERE the bores sit, not just how big they are. Two blocks drilled in
    # different places have the same outside size, the same radii and the same
    # VOLUME -- the holes are identical, only their positions differ -- so
    # without this they merge and every copy inherits the first one's hole
    # pattern. Shape 9 was the case: solids 200 and 222 would have been built
    # with solid 43's bores, 203 mm from where they belong. Their principal
    # moments of inertia differ by 7%, which is what proves they are genuinely
    # different parts rather than one part measured from the other end.
    #
    # Measured in the solid's OWN frame, as unsigned distances from its centre,
    # so the key stays free of position and orientation: a brick turned
    # end-for-end still matches itself, while a brick drilled elsewhere does
    # not.
    bores = tuple(sorted(bore_offsets(solid)))
    return (kind, sig, dims, radii, bores, round(solid.Volume, 1))


def bore_offsets(solid):
    """Each bore's distance from the solid's centre along its own axes.

    Unsigned, because the sign depends on which end of the block the frame
    happens to point at, and that is genuinely arbitrary -- it is the only part
    of the offset a rotation can change for an axis-aligned part.

    Rounded to the NEAREST MILLIMETRE, not to a fraction of one. The de-tilt
    leaves a few hundredths of a millimetre of scatter on every measurement,
    so a finer quantisation splits identical parts: solids 6 and 108 read 9.53
    against 9.52 and were filed as different shapes, when a boolean test shows
    them overlapping to 100.000000%. A whole millimetre is far below the
    spacing of genuinely different hole patterns -- shape 9's two patterns are
    203 mm apart -- and far above the noise.
    """
    axes = orthogonal_frame(solid)
    cyls = cylinders(solid)
    if axes is None or not cyls:
        return []
    c = solid.BoundBox.Center
    out = []
    for f in cyls:
        cc = f.Surface.Center
        d = FreeCAD.Vector(cc.x - c.x, cc.y - c.y, cc.z - c.z)
        out.append(tuple(sorted(round(abs(d.dot(u))) for u in axes))
                   + (round(f.Surface.Radius, 1),))
    return out


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

# The CAD label is the best name we have, so names are built from it rather
# than from the colour. Two fixes are applied on the way:
#
#   * a label naming the wrong material. The colour-derived names were worse
#     than unhelpful -- shape 5 was "SteelBrick1x6x24PV" for a BP poly edge,
#     shape 19 "SteelRwall19PV" for lead, shape 25 "SteelTube25PV" for copper.
#   * "Poly Lead Block" reads as a lead part. It is polyethylene sitting in the
#     lead, so it is renamed to avoid exactly that confusion.
# The description part of the name, keyed by CAD label. Everything here is a
# deliberate edit to what the label says -- fixing a typo, dropping a word that
# carries no meaning, or naming the part for what it is. Labels not listed are
# used as written.
LABEL_RENAME = {
    # Poly, but sitting in the lead: the label reads as a lead part.
    "Poly Lead Block": "Insert Block",
    "Poly Block Reduced": "Block",
    # "SCC" is a typo for SSC, and these are the steel blocks around it.
    "Steel SCC Blocks": "SSC Blocks",
    # "S sq" means nothing in a volume name.
    "Poly Top S sq walls": "Top Square Wall",
    "Top Al sq walls": "Top Square Wall",
    # The two 2-inch pipes, named by what distinguishes them: length.
    "LaBr Cu": "Long Pipe",
    "Cu Insert Front": "Short Pipe",
    # "Changed" and "Machining" describe a CAD operation, not a part.
    "Shelf Changed": "Shelf",
    "BASE PLATE MACHINING": "Base Plate",
    # Size codes read better spelled out: 2x4x8 inch, not 248.
    "248 Lead Brick": "Brick 2x4x8",
    "2416 Lead Brick": "Brick 2x4x16",
    # Dropping the redundant material word would leave these with no noun at
    # all -- "SSC Poly" becomes bare "SSC", "Front Side Small Poly" becomes
    # "Front Side Small". Give them one.
    "SSC Poly": "SSC Panel",
    "Front Side Small Poly": "Front Side Small Panel",
    # The CAD labels use both "Bot" and "Bottom"; spell it out everywhere.
    # "Cu Bottom Cut" needs no entry -- it already says Bottom, and the "Cu"
    # is dropped as redundant with the material prefix.
    "Small SSC Bot": "Small SSC Bottom",
    "Bot Square Poly": "Bottom Square",
}

# Words that merely repeat the material prefix, dropped so the name does not
# read "BPPolyFrontOuter". Matched case-insensitively as whole words.
_REDUNDANT = ("lead", "poly", "cu", "al", "aluminum", "steel", "tungsten",
              "bp")

# Material words, for labels that carry none ("Shelf Changed", "Small SSC Top")
# and for labels that start with a digit ("248 Lead Brick"), which cannot lead
# a C++ identifier.
MATERIAL_NAME = {
    "G4_Pb": "Lead",
    "BP": "BP",                 # borated polyethylene -- NOT "Poly"
    "Polyethylene": "Poly",     # the two genuinely plain-poly parts
    "CollCu": "Cu",
    "Al6061": "Aluminum",
    "MildSteel": "Steel",
    "G4_W_Hayman": "Tungsten",
}

def label_name(raw, material):
    """Volume name from the CAD label: 'F10256760--_1-248 Lead Brick_001'.

    One scheme, applied to every part: {Material}{Description}PV, material
    always first. The first pass took whatever order the label happened to use
    and was incoherent as a result -- PolyFrontOuterPV, SmallTriangleLeadPV and
    TopAlSqWallsPV put the material at the front, the end and the middle.

    The material comes from the part label (see solid_materials.json), so BP
    parts say BP rather than Poly: of the 22 poly parts only two are plain
    polyethylene, and calling the other twenty "Poly" hid that.

    A word in the description that merely repeats the material is dropped, so
    "Poly Front Outer" in BP becomes BPFrontOuter, not BPPolyFrontOuter.

    No shape id is appended: most labels are already unique, and a numeric
    suffix runs into any trailing size code. The caller disambiguates the few
    labels genuinely shared by two shapes.

    Returns None when there is no usable label, so the caller can fall back.
    """
    if not raw:
        return None
    s = re.sub(r"^F\d+-*", "", raw)          # part number
    s = re.sub(r"^_\d+-", "", s)             # revision prefix
    s = re.sub(r"_\d+$", "", s)              # instance suffix
    s = LABEL_RENAME.get(s.strip(), s.strip())
    if not s:
        return None

    word = MATERIAL_NAME.get(material, "")

    # Title-case each word, so an all-caps label ("BASE PLATE MACHINING") does
    # not survive as an unreadable run of capitals. Words that are already
    # mixed case are left alone -- "SSC", "LaBr" and "Cu" are meaningful as
    # written and lower-casing them would lose that.
    parts = []
    for w in re.split(r"[^A-Za-z0-9]+", s):
        if not w:
            continue
        parts.append(w.capitalize() if w.isupper() and len(w) > 3 else
                     w[:1].upper() + w[1:])

    # Drop words that only restate the material, then put the material first.
    # A size code keeps its place at the end of the description ("2x4x8"),
    # which is also what keeps the name a valid identifier.
    parts = [w for w in parts if w.lower() not in _REDUNDANT]
    if word:
        parts.insert(0, word)

    return "".join(parts) + "PV"


def suggest_name(kind, p, colour, seq, count, spread):
    """A starting point for the volume name, following the Offline scheme.

    Used only when the solid carries no CAD label. Offline names shielding as
    {Material}{Side}wall{N}PV -- CopperLwallPV, LeadTwall1PV, BPRwall2PV. That
    only reads well for a panel that sits on one side of the house. Most of
    this model is stacked brick: 146 copies of one 2x4x8 inch block scattered
    over more than a metre in every direction, where a side label would be
    actively misleading. So repeated, widely spread shapes are named for what
    they are, and only shapes that stay put get a side.
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
        key = v2_key(i)
        if key is None:
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

    # Name every shape first, so collisions can be seen and broken.
    #
    # Four labels are each used by two DIFFERENT shapes -- two "SSC Poly"
    # panels of different size, two "Cu Bottom Cut", and so on. Appending the
    # shape id to every name to guard against that made the other 38 worse
    # ("LeadBrick2416" + "0" reads as "LeadBrick24160"), so the suffix goes
    # only where it is earned, and as a letter to keep it out of the digits.
    shape_names = {}
    for key, g in groups.items():
        first = g["members"][0]
        nm = label_name(labels.get(first), SOLID_MATERIALS.get(str(first)))
        if not nm:
            pts = [placements[i][2]["pos"] for i in g["members"]]
            spread = max(max(q[a] for q in pts) - min(q[a] for q in pts)
                         for a in range(3)) if len(pts) > 1 else 0.0
            nm = suggest_name(g["kind"], g["p"], per_solid.get(first, ""),
                              order[key], len(g["members"]), spread)
        shape_names[key] = nm
    clash = collections.Counter(shape_names.values())
    seen = collections.Counter()
    for key in list(shape_names):
        nm = shape_names[key]
        if clash[nm] > 1:
            suffix = chr(ord("A") + seen[nm])
            seen[nm] += 1
            shape_names[key] = "%s%sPV" % (nm[:-2], suffix)

    shape_rows = []
    for key, g in groups.items():
        n = order[key]
        p = g["p"]
        dims = (p.get("dx", 0.0), p.get("dy", 0.0), p.get("dz", 0.0))

        # Iteration 2's dimensions, for a group that came from the v2 map.
        #
        # Identity from v2, geometry from v3, SIZE from v2: iteration 3 shrank
        # most parts within manufacturing tolerance, so its measured sizes are
        # 24 subtly different lead bricks where there is physically one. Taking
        # the size from v2 restores the single definition.
        #
        # Only dx/dy/dz are overridden. The prism and tube fields are left on
        # v3 deliberately: the cap outline written to stm_prisms.csv comes from
        # the v3 solid, so a v2 cap_area/sweep_len here would contradict the
        # outline file. Solid 220 is the case that matters -- v2 identity, but
        # v3's 8-face geometry.
        v2 = None
        for _m in g["members"]:
            v2 = v2_dims(_m)
            if v2:
                break
        if v2 and v2.get("dx") not in (None, ""):
            try:
                dims = (float(v2["dx"]), float(v2["dy"]), float(v2["dz"]))
            except (TypeError, ValueError):
                pass

        # Colour of the group: the members agree, so take the commonest and
        # fall back to the dimension lookup only if none of them matched.
        votes = collections.Counter(per_solid[i] for i in g["members"]
                                    if i in per_solid)
        if votes:
            colour = votes.most_common(1)[0][0]
            how = "part" if len(votes) == 1 else "part(mixed)"
        else:
            colour, how = colour_from_lut(dims, lut)

        # Name from the CAD label when there is one: it says what the part
        # actually is, where the colour-derived fallback could only guess from
        # size and position. suggest_name() stays for solids with no label.
        name = shape_names[key]

        # dx/dy/dz are the SOLID's dimensions for a box, a bored box or a tube.
        # For a prism they would be the bounding envelope, which is not the
        # shape: the solid is a swept polygon that fills only part of that box,
        # and reading the envelope as a G4Box overstates it (shape 19 is a
        # wedge of 2903 mm^2 cap area inside a 76 x 152 x 76 envelope, so a box
        # would be twice the material). The envelope stays available from the
        # placements and from stm_prisms.csv, so nothing is lost by leaving
        # these blank; what is gained is that the columns cannot be mistaken
        # for a buildable size.
        is_prism = "prism_outline" in p
        row = {
            "shape_id": n,
            "type": g["kind"],
            "count": len(g["members"]),
            "name": name,
            "material": colour,
            "material_hint": COLOUR_HINT.get(colour, ""),
            "material_src": how,
            "dx": "" if is_prism else round(dims[0], 3),
            "dy": "" if is_prism else round(dims[1], 3),
            "dz": "" if is_prism else round(dims[2], 3),
            "rmin": round(p["rmin"], 3) if "rmin" in p else "",
            "rmax": round(p["rmax"], 3) if "rmax" in p else "",
            "nbores": len(p.get("bores", [])) or "",
            "rotY45": "yes" if g["oblique"] else "",
            "in_x": "" if is_prism else round(dims[0] / 25.4, 3),
            "in_y": "" if is_prism else round(dims[1] / 25.4, 3),
            "in_z": "" if is_prism else round(dims[2] / 25.4, 3),
            # The prism's own defining numbers, so this row is self-contained:
            # a cap of this area swept this far, with the outline in
            # stm_prisms.csv. Blank for every other kind.
            "cap_area": round(p["prism_area"], 3) if is_prism else "",
            "sweep_len": round(p["prism_len"], 3) if is_prism else "",
            "cap_verts": len(p["prism_outline"]) if is_prism else "",
        }
        # A PRISM that failed the extrusion test keeps its envelope, since that
        # is genuinely all that is known about it. Nothing in this model does,
        # but a future export might, and silently blanking the only dimensions
        # available would hide the part rather than describe it.
        if g["kind"] in ("PRISM", "OTHER") and not is_prism:
            row["dx"] = round(dims[0], 3)
            row["dy"] = round(dims[1], 3)
            row["dz"] = round(dims[2], 3)
            row["in_x"] = round(dims[0] / 25.4, 3)
            row["in_y"] = round(dims[1] / 25.4, 3)
            row["in_z"] = round(dims[2] / 25.4, 3)
        shape_rows.append(row)

    # One row per solid: where each block goes, in Mu2e coordinates relative
    # to _STMShieldingRef, plus its bores if it has any.
    place_rows = []
    for i, key, p, s in placements:
        x, y, z = p["pos"]
        # Orientation of this copy, as a 3x3 rotation taking the shape's own
        # frame to Mu2e axes. Without it the CSVs are not self-contained: two
        # copies of one shape can sit at right angles and differ in no recorded
        # column, so a consumer cannot place either. Shape 1 has 13 distinct
        # orientations across its 146 copies.
        # The rotation the CLASSIFIER computed, alongside the canonical
        # dimensions it reported. Re-deriving it here from the solid's face
        # normals -- which is what solid_orientation() did -- produced a
        # rotation relative to the solid's own axis ORDER while dx/dy/dz were
        # ordered canonically, so composing the two double-counted the
        # orientation. Reading it from the same call that fixed the frame is
        # what keeps them consistent.
        rot = p.get("rot")
        place_rows.append({
            "solid_id": i,
            "shape_id": order[key],
            "part": labels.get(i, ""),
            "x": round(x, 3),
            "y": round(y, 3),
            "z": round(z, 3),
            # Row-major 3x3: world_axis = sum_j r{row}{col} * own_axis.
            # "" for a solid with no axis-aligned orientation (a wedge, or a
            # 45 deg part) -- reported rather than guessed at.
            "r11": rot[0][0] if rot else "", "r12": rot[0][1] if rot else "",
            "r13": rot[0][2] if rot else "", "r21": rot[1][0] if rot else "",
            "r22": rot[1][1] if rot else "", "r23": rot[1][2] if rot else "",
            "r31": rot[2][0] if rot else "", "r32": rot[2][1] if rot else "",
            "r33": rot[2][2] if rot else "",
            "orientation": orientation_label(rot),
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
                # Position RELATIVE to the block centre: this is what a
                # G4SubtractionSolid needs, and unlike an absolute position it
                # is valid for every placement of the shape, not just the one
                # the bore happened to be measured from.
                "dx": round(b["rel"][0], 3),
                "dy": round(b["rel"][1], 3),
                "dz": round(b["rel"][2], 3),
                "depth": round(b["depth"], 3),
                # Absolute Mu2e position too, for cross-checking against the
                # CAD. Derived, not authoritative: dx/dy/dz are the definition.
                "x": round(b["pos"][0], 3),
                "y": round(b["pos"][1], 3),
                "z": round(b["pos"][2], 3),
                "ax": round(b["axis"][0], 4),
                "ay": round(b["axis"][1], 4),
                "az": round(b["axis"][2], 4),
            })

    # One row per vertex of every prism's cap outline, in order. Together with
    # the axis and length this is a complete G4ExtrudedSolid: the CSVs no longer
    # stop at a bounding envelope for these parts.
    prism_rows = []
    for key, g in groups.items():
        p = g["p"]
        if "prism_outline" not in p:
            continue
        n = order[key]
        ax = p["prism_axis"]
        u = p["prism_u"]
        w = p["prism_w"]
        for seq, (a, b2) in enumerate(p["prism_outline"]):
            prism_rows.append({
                "shape_id": n,
                "seq": seq,
                # Cap outline in the cap's own plane, mm. Sweep this polygon
                # along prism_axis for prism_len to rebuild the solid.
                "u": a,
                "v": b2,
                "len": round(p["prism_len"], 3),
                "area": round(p["prism_area"], 3),
                # The cap plane's basis and the sweep direction, in Mu2e
                # coordinates, so the 2D outline can be placed in 3D.
                # 10 decimal places, not the 4 used elsewhere. These are unit
                # vectors multiplied by an arm of up to ~250 mm, so a rounding
                # of 1e-4 rad displaces a rebuilt vertex by ~0.025 mm. At 1e-10
                # the round-trip is limited by the outline and position columns
                # instead, i.e. by real geometry rather than by formatting.
                "axis_x": round(ax[0], 10),
                "axis_y": round(ax[1], 10),
                "axis_z": round(ax[2], 10),
                "u_x": round(u[0], 10),
                "u_y": round(u[1], 10),
                "u_z": round(u[2], 10),
                "w_x": round(w[0], 10),
                "w_y": round(w[1], 10),
                "w_z": round(w[2], 10),
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
    if prism_rows:
        write(os.path.join(OUTDIR, "stm_prisms.csv"), prism_rows)

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
