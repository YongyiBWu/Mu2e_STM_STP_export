"""
Plot every shape to its own PNG, with the Geant4 placement point marked.

    "C:/Program Files/FreeCAD 1.1/bin/python.exe" scripts/plot_all_components.py
    "C:/Program Files/FreeCAD 1.1/bin/python.exe" scripts/plot_all_components.py --outdir plots

Writes one PNG per shape into plots/ (default), named
shape_<id>_<TYPE>_<name>.png, plus a copy of plot_tilt_state.py's figure so the
directory holds the whole picture of the model in one place.

What "location" means here, per shape kind
    BOX, BOX_HOLE   the solid's centre. A G4Box is centred on its origin, so
                    the placement in stm_placements.csv IS the marked point.
    TUBE            the tube's centre, likewise -- G4Tubs is centred on its
                    own axis midpoint.
    PRISM           the (0,0) of the cap plane, at the sweep half-length.

    That last one is the Geant4 convention for G4ExtrudedSolid: the solid is
    built from a 2D polygon swept +/-len/2 about z in its own frame, and its
    origin is wherever (0,0) falls in that polygon -- NOT the polygon's
    centroid. stm_prisms.csv stores the outline relative to the solid's centre,
    so (0,0) is that centre, and the placement from stm_placements.csv is
    exactly the point to hand G4PVPlacement.

    The distinction matters: for shape 13 the cap's area centroid sits 72 mm
    from (0,0). Marking the centroid would put the label somewhere Geant4 will
    not place anything.

Every plot marks the placement point with a crosshair and prints its Mu2e
coordinates, so the picture and the CSV can be checked against each other.
"""

import argparse
import collections
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
OUTDIR = os.path.join(ROOT, "output")

sys.path.insert(0, HERE)
import show_shape as S                                    # noqa: E402


def safe(name):
    """A filename fragment with nothing a shell or filesystem will object to."""
    keep = [c if (c.isalnum() or c in "-_") else "_" for c in (name or "")]
    return "".join(keep).strip("_") or "unnamed"


def placement_point(shape, placements):
    """Where Geant4 places this shape, in Mu2e mm, and what that point is.

    For every kind this is the solid's own origin, which is what the
    placements CSV records. The value of spelling it out per kind is that the
    caller can print WHY the point is what it is -- a prism's origin is the
    cap plane's (0,0), which is not its centroid and not obviously its centre.
    """
    if not placements:
        return None, "no placement row"
    p = placements[0]
    pt = (S.num(p, "x", 0.0), S.num(p, "y", 0.0), S.num(p, "z", 0.0))
    kind = shape["type"]
    if kind == "PRISM":
        why = "cap-plane (0,0) at sweep half-length (G4ExtrudedSolid origin)"
    elif kind == "TUBE":
        why = "tube centre on its axis (G4Tubs origin)"
    else:
        why = "block centre (G4Box origin)"
    return pt, why


def placement_lines(placements, limit=8):
    """Every copy's position, not just the first.

    A shape with 146 copies cannot list them all on a figure, so the tail is
    summarised rather than dropped -- the full list lives in
    output/placements_by_shape.csv, and the footnote says so.
    """
    n = len(placements)
    if n == 0:
        return "no placement rows"
    head = ["placed at Mu2e (x, y, z) mm, with each copy's rotation"]
    head.append("from the canonical frame  --  %d placement(s):" % n)
    for p in placements[:limit]:
        head.append("   solid %-4s %-16s (%9.3f, %9.3f, %9.3f)  %s"
                    % (p["solid_id"], p.get("part", "?"),
                       S.num(p, "x", 0.0), S.num(p, "y", 0.0),
                       S.num(p, "z", 0.0),
                       (p.get("orientation") or "(no rotation)")))
    if n > limit:
        head.append("   ... %d more; full list in output/placements_by_shape.csv"
                    % (n - limit))
    return "\n".join(head)


def model_boxes():
    """Every solid's Mu2e-frame bounding box, for the context panel.

    Built from stm_placements.csv + stm_shapes.csv rather than
    output/tilt_state.json: that file holds CAD-frame centres on purpose (it
    describes the input geometry), so using it here would draw every part
    hundreds of mm from where it actually sits.

    A prism has no dx/dy/dz by design, so its extent comes from the sweep
    length and the cap's bounding size in the outline plane.
    """
    shapes = {r["shape_id"]: r for r in S.load("stm_shapes.csv")}
    prisms = collections.defaultdict(list)
    for r in S.load("stm_prisms.csv", required=False):
        prisms[r["shape_id"]].append(r)

    out = []
    for p in S.load("stm_placements.csv"):
        sh = shapes.get(p["shape_id"], {})
        dx = S.num(sh, "dx")
        if dx is None:
            rows = prisms.get(p["shape_id"], [])
            if not rows:
                continue
            # Build the prism's world extent from its ACTUAL geometry: the cap
            # spans u and v along the stored basis vectors, and the sweep runs
            # `len` along the axis. Projecting those three onto the world axes
            # gives the real bounding box.
            #
            # The previous version took one span = max(u_extent, v_extent),
            # used it for BOTH dx and dy, and put len on dz unconditionally.
            # That squared off a rectangular cap and pinned the sweep to z
            # whatever direction it really ran -- shape 13 has a 490 x 815 cap
            # swept 25.4 along X, and was drawn as 815 x 815 x 25.4, the wrong
            # shape in the wrong orientation.
            r0 = rows[0]
            us = [S.num(r, "u", 0.0) for r in rows]
            vs = [S.num(r, "v", 0.0) for r in rows]
            du = max(us) - min(us)
            dv = max(vs) - min(vs)
            ln = S.num(r0, "len", 0.0)
            uvec = (S.num(r0, "u_x", 1.0), S.num(r0, "u_y", 0.0),
                    S.num(r0, "u_z", 0.0))
            wvec = (S.num(r0, "w_x", 0.0), S.num(r0, "w_y", 1.0),
                    S.num(r0, "w_z", 0.0))
            avec = (S.num(r0, "axis_x", 0.0), S.num(r0, "axis_y", 0.0),
                    S.num(r0, "axis_z", 1.0))
            # Extent along each world axis: each contributing vector projects
            # onto it, and for an axis-aligned prism exactly one dominates.
            dims = tuple(abs(uvec[k]) * du + abs(wvec[k]) * dv
                         + abs(avec[k]) * ln for k in range(3))
        else:
            # Rotate the CANONICAL dimensions into world axes for this copy.
            #
            # stm_shapes.csv stores dx/dy/dz in the shape's canonical frame --
            # longest first for a box -- so using them directly as world
            # extents draws every rotated copy with the wrong profile. Solid 0
            # is 406.4 x 101.6 x 50.8 canonically but occupies
            # 406.4 x 50.8 x 101.6 in the model, and solid 7 of the same shape
            # occupies 101.6 x 406.4 x 50.8. Before the canonical rework the
            # stored dims happened to be world-ordered and this was right by
            # accident.
            #
            # |R| applied to the extents gives the rotated bounding box: each
            # world axis takes whichever canonical extent the rotation maps
            # onto it. A copy with no recorded rotation keeps its canonical
            # dims, which is the best available answer.
            cd = (dx, S.num(sh, "dy", 0.0), S.num(sh, "dz", 0.0))
            rot = rotation_of(p)
            if rot is None:
                dims = cd
            else:
                dims = tuple(sum(abs(rot[k][j]) * cd[j] for j in range(3))
                             for k in range(3))
        out.append({"solid_id": p["solid_id"], "shape_id": p["shape_id"],
                    "x": S.num(p, "x", 0.0), "y": S.num(p, "y", 0.0),
                    "z": S.num(p, "z", 0.0),
                    "dx": dims[0], "dy": dims[1], "dz": dims[2]})
    return out


def add_context(fig, placements, boxes, chosen_id=None):
    """Inset showing where this shape's copies sit in the whole model.

    Without it a plot says what a part looks like but not where it belongs,
    which is the first thing anyone asks of an unfamiliar part number. The
    whole house is drawn faint; this shape's copies are drawn solid on top.
    """
    if not boxes:
        return
    mine = {p["solid_id"] for p in placements}
    # Draw into the axes show_shape.draw() reserved at the top of the left
    # column. Creating a floating add_axes rectangle here instead -- which is
    # what this did before the fixed layout -- puts the inset at a position
    # that depends on the canvas and can land on top of a title.
    ax = getattr(fig, "_stm_context", None)
    if ax is None:
        return
    ax.set_facecolor("none")

    def edges(b):
        hx, hy, hz = b["dx"] / 2.0, b["dy"] / 2.0, b["dz"] / 2.0
        x, y, z = b["x"], b["y"], b["z"]
        c = [(x - hx, y - hy, z - hz), (x + hx, y - hy, z - hz),
             (x + hx, y + hy, z - hz), (x - hx, y + hy, z - hz),
             (x - hx, y - hy, z + hz), (x + hx, y - hy, z + hz),
             (x + hx, y + hy, z + hz), (x - hx, y + hy, z + hz)]
        for a, b2 in [(0, 1), (1, 2), (2, 3), (3, 0), (4, 5), (5, 6),
                      (6, 7), (7, 4), (0, 4), (1, 5), (2, 6), (3, 7)]:
            # Mu2e (x, y, z) -> matplotlib (horizontal, depth, vertical), so
            # the inset shares the +y-up convention of the panel beside it.
            pa, pb = S.mu2e_axes(*c[a]), S.mu2e_axes(*c[b2])
            yield ([pa[0], pb[0]], [pa[1], pb[1]], [pa[2], pb[2]])

    for b in boxes:
        hit = b["solid_id"] in mine
        if hit:
            continue
        for xs, ys, zs in edges(b):
            ax.plot(xs, ys, zs, color="#b0b7bd", lw=0.25, alpha=0.45)
    # With 146 copies the solid outlines merge into a red mass that hides where
    # they actually are, so thin the stroke as the count grows.
    hits = [b for b in boxes if b["solid_id"] in mine]
    lw = 1.1 if len(hits) <= 8 else (0.7 if len(hits) <= 40 else 0.45)
    for b in hits:
        # The one copy drawn in the "as placed" pad gets that pad's colour, so
        # a shape with 20 copies still says WHICH copy the rotated view shows.
        is_chosen = chosen_id is not None and b["solid_id"] == chosen_id
        ax.plot([], [], [])
        for xs, ys, zs in edges(b):
            ax.plot(xs, ys, zs,
                    color=PLACED_COLOUR if is_chosen else "#c0392b",
                    lw=lw * 2.0 if is_chosen else lw,
                    alpha=1.0 if is_chosen else 0.95,
                    zorder=6 if is_chosen else 4)

    lim = {}
    for k in "xyz":
        lo = min(b[k] - b["d" + k] / 2.0 for b in boxes)
        hi = max(b[k] + b["d" + k] / 2.0 for b in boxes)
        lim[k] = (lo, hi)
    # Axis slots hold Mu2e x, z, y in that order (see mu2e_axes).
    ax.set_xlim(*lim["x"])
    ax.set_ylim(*lim["z"])
    ax.set_zlim(*lim["y"])
    try:
        ax.set_box_aspect([lim["x"][1] - lim["x"][0],
                           lim["z"][1] - lim["z"][0],
                           lim["y"][1] - lim["y"][0]])
    except (AttributeError, ValueError):
        pass
    # Same Mu2e orientation AND the same raised viewpoint as the part panel,
    # so the inset is neither a mirror image nor a different angle from the
    # plot it sits beside.
    ax.view_init(elev=32, azim=-68)
    ax.set_xlim(ax.get_xlim()[::-1])
    for a in (ax.xaxis, ax.yaxis, ax.zaxis):
        a.set_ticklabels([])
    ax.grid(False)
    ax.set_title("location in model (red = this shape)\n"
                 "+y up, -x right, +z into page", fontsize=6.5,
                 color="#444444")


# The copy shown in the "as placed" pad, and highlighted in the context view.
# One colour, used in both, so the two panels are unambiguously about the same
# solid -- otherwise "which of the 20 copies is this?" has no answer on a page
# showing a shape with many placements.
PLACED_COLOUR = "#1a7f37"


def rotation_of(place):
    """This copy's canonical -> Mu2e rotation as a 3x3, or None."""
    if not (place.get("r11") or "").strip():
        return None
    return [[S.num(place, "r%d%d" % (i + 1, j + 1), 0.0) for j in range(3)]
            for i in range(3)]


def _bore_faces(bores, shape, S_):
    """Cylinder faces for every bore, in the shape's canonical frame.

    The rotated pad showed the outer block only, which is misleading for a
    bored part: the pad exists to show what the solid actually looks like once
    placed, and a block with four holes through it does not look like a block.

    Offsets and axes in stm_bores.csv are already canonical, so they need no
    conversion here -- they compose with the box the pad draws.
    """
    out = []
    dx = S_.num(shape, "dx", 0.0)
    dy = S_.num(shape, "dy", 0.0)
    dz = S_.num(shape, "dz", 0.0)
    for b in bores:
        r = S_.num(b, "r", 0.0)
        c = (S_.num(b, "dx", 0.0), S_.num(b, "dy", 0.0), S_.num(b, "dz", 0.0))
        ax = (S_.num(b, "ax", 0.0), S_.num(b, "ay", 0.0), S_.num(b, "az", 1.0))
        depth = S_.num(b, "depth", 0.0)
        if depth <= 0:
            depth = abs(ax[0]) * dx + abs(ax[1]) * dy + abs(ax[2]) * dz
        out.extend(S_.cylinder_faces(c, ax, r, depth))
    return out


def add_placed(fig, shape, placements, geom, bores):
    """Draw the shape AFTER rotation, for one representative copy.

    The main panel draws the canonical frame, because that is the frame the
    dx/dy/dz describe. That is correct but it hides the orientation: two copies
    of one shape can sit at right angles and look identical on this page. This
    pad closes the gap by showing what the part actually looks like once its
    rotation is applied, and names the copy it belongs to.

    The copy chosen is the first one that HAS a rotation, so a shape whose
    first placement is oblique still gets a meaningful pad.
    """
    ax = getattr(fig, "_stm_placed", None)
    if ax is None:
        return
    import matplotlib.pyplot as plt                       # noqa: F401
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection

    chosen, rot = None, None
    for p in placements:
        r = rotation_of(p)
        if r is not None:
            chosen, rot = p, r
            break
    if chosen is None:
        chosen = placements[0] if placements else None

    # Faces in the canonical frame, from the same helpers the main panel uses.
    bore_faces = []
    if geom is not None:
        faces = S.prism_faces(geom)
    else:
        dx = S.num(shape, "dx", 0.0)
        dy = S.num(shape, "dy", 0.0)
        dz = S.num(shape, "dz", 0.0)
        if shape["type"] == "TUBE":
            rmax = S.num(shape, "rmax", max(dx, dy) / 2.0)
            faces = S.cylinder_faces((0.0, 0.0, 0.0), (0, 0, 1), rmax,
                                     S.num(shape, "dz", 0.0))
        else:
            faces = S.box_faces(0.0, 0.0, 0.0, dx, dy, dz)
            # A bored block is not a block. Carry the holes into this view too,
            # or the pad claims the part is solid on a page whose section
            # panels are entirely about its bores.
            bore_faces = _bore_faces(bores, shape, S)

    # Every kind now arrives in its LOCAL frame -- prism_faces() emits the cap
    # in x,y swept along z, matching box_faces() and cylinder_faces() -- so the
    # rotation applies uniformly. Suppressing it for prisms was only correct
    # while prism_faces() returned world-placed geometry; with that removed,
    # skipping the rotation here would show the local frame in both panels and
    # leave the placement invisible.
    already_placed = False

    # S.* helpers already emit matplotlib axis order (x, z, y); undo that to
    # get Mu2e order, rotate, then re-apply. Rotating the display order would
    # silently transpose two axes.
    def rotated(pt):
        x, z, y = pt                       # as emitted by mu2e_axes
        if rot is None or already_placed:
            return S.mu2e_axes(x, y, z)
        v = (x, y, z)
        w = tuple(sum(rot[k][j] * v[j] for j in range(3)) for k in range(3))
        return S.mu2e_axes(*w)

    rfaces = [[rotated(p) for p in f] for f in faces]
    ax.add_collection3d(Poly3DCollection(
        rfaces, alpha=0.35 if bore_faces else 0.45, facecolor=PLACED_COLOUR,
        edgecolor="#14532d", linewidth=0.8))

    # Bores on top of the block, in the same red the section panels use, so a
    # reader moving between the two panels is looking at the same holes. The
    # block goes slightly more transparent above when there are bores, or they
    # are hidden inside it.
    rbores = [[rotated(p) for p in f] for f in bore_faces]
    if rbores:
        ax.add_collection3d(Poly3DCollection(
            rbores, alpha=0.95, facecolor="#c0392b", edgecolor="#7b241c",
            linewidth=0.3))

    pts = [p for f in rfaces for p in f] + [p for f in rbores for p in f]
    for setter, k in ((ax.set_xlim, 0), (ax.set_ylim, 1), (ax.set_zlim, 2)):
        lo = min(p[k] for p in pts)
        hi = max(p[k] for p in pts)
        pad = max((hi - lo) * 0.08, 1e-6)
        setter(lo - pad, hi + pad)
    try:
        ax.set_box_aspect([max(p[k] for p in pts) - min(p[k] for p in pts)
                           or 1.0 for k in range(3)])
    except (AttributeError, ValueError):
        pass
    ax.view_init(elev=32, azim=-68)
    ax.set_xlim(ax.get_xlim()[::-1])
    for a in (ax.xaxis, ax.yaxis, ax.zaxis):
        a.set_ticklabels([])
    ax.grid(False)
    if chosen is not None:
        ax.set_title("as placed: solid %s\n%s"
                     % (chosen["solid_id"],
                        chosen.get("orientation") or "no rotation recorded"),
                     fontsize=7, color="#14532d")
    return chosen


def _prism_rows_for(shape):
    """The stm_prisms.csv rows for this shape, or [] if it is not a prism.

    Loaded here rather than threaded through every caller: the excerpt is the
    only consumer, and a prism's outline is small.
    """
    rows = [r for r in S.load("stm_prisms.csv", required=False)
            if r["shape_id"] == shape["shape_id"]]
    return sorted(rows, key=lambda r: int(r["seq"]))


def g4_excerpt(shape, placements, chosen):
    """The Geant4 construction for ONE copy, for the text pad.

    The same code scripts/emit_g4.py writes to
    output/constructSTM_generated.cc, narrowed to the copy the "as placed" pad
    shows. Putting it on the page means the plot, the rotation column and the
    code a consumer will actually write can be read against each other without
    opening three files.

    Kept deliberately short -- the solid plus one placement. The full file has
    every copy, the bore subtractions and the material calls.
    """
    try:
        sys.path.insert(0, HERE)
        import emit_g4
    except Exception:
        return ""
    if chosen is None:
        chosen = placements[0] if placements else None
    if chosen is None:
        return ""

    var = "s%s" % shape["shape_id"]
    lines = ["G4 (one copy; full file: output/constructSTM_generated.cc)"]

    # For a prism the polygon IS the shape, so show every vertex rather than
    # trimming to the first few lines. For a box the first statement carries
    # the dimensions and the rest is bore subtraction, which is summarised.
    prism_rows = _prism_rows_for(shape)
    if prism_rows:
        solid = emit_g4.solid_code(shape, [], prism_rows, var)
        for ln in solid.splitlines():
            s = ln.strip()
            if s.startswith("//"):
                continue                   # the basis commentary, not needed
            lines.append("  " + s)
    else:
        solid = emit_g4.solid_code(shape, [], [], var)
        # Only the first statement of the solid: bored blocks add a G4Tubs and
        # a subtraction per bore, which does not fit and is not the point here.
        for ln in solid.splitlines()[:3]:
            lines.append("  " + ln.strip())
    # Say so explicitly. Without this the excerpt for a BOX_HOLE reads as a
    # plain G4Box on a page whose other two panels are entirely about its
    # holes -- an omission that looks like a contradiction.
    nb = len([b for b in (shape.get("nbores") or "") if b.isdigit()])
    if (shape.get("nbores") or "").strip():
        # Say clearly that the FILE has them and only this excerpt is short.
        # "omitted here" read as though the emitter dropped the bores.
        lines.append("  // %s G4Tubs subtraction(s) follow in the generated"
                     % shape["nbores"])
        lines.append("  // file -- trimmed from this excerpt only, for space")
    rot = emit_g4.rotation_code(chosen, "%s_rot" % var)
    for ln in rot.splitlines():
        lines.append("  " + ln.strip())
    lines.append("  new G4PVPlacement(%s_rot," % var)
    lines.append("      G4ThreeVector(%.3f, %.3f, %.3f)*CLHEP::mm,"
                 % (S.num(chosen, "x", 0.0), S.num(chosen, "y", 0.0),
                    S.num(chosen, "z", 0.0)))
    lines.append("      %sLV, \"%s_%s\", motherLV, false, %s);"
                 % (var, emit_g4.ident(shape.get("name")),
                    chosen["solid_id"], chosen["solid_id"]))
    return "\n".join(lines)


def mark_origin(ax, lim):
    """Crosshair at the solid's origin, which is where Geant4 places it.

    Drawn as three axis-spanning lines rather than a dot: on a 900 mm slab a
    dot at the centre is invisible, and the lines also show the frame the
    dimensions are measured in.
    """
    ax.plot([lim[0][0], lim[0][1]], [0, 0], [0, 0],
            color="#c0392b", lw=1.0, alpha=0.8, zorder=10)
    ax.plot([0, 0], [lim[1][0], lim[1][1]], [0, 0],
            color="#c0392b", lw=1.0, alpha=0.8, zorder=10)
    ax.plot([0, 0], [0, 0], [lim[2][0], lim[2][1]],
            color="#c0392b", lw=1.0, alpha=0.8, zorder=10)
    ax.scatter([0], [0], [0], color="#c0392b", s=28, depthshade=False,
               zorder=11)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default=os.path.join(ROOT, "plots"))
    ap.add_argument("--dpi", type=int, default=120)
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    shapes = S.load("stm_shapes.csv")
    places = S.load("stm_placements.csv")
    bores_all = S.load("stm_bores.csv")
    prisms_all = S.load("stm_prisms.csv", required=False)
    all_boxes = model_boxes()

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    print("writing %d shape plots to %s" % (len(shapes), args.outdir))
    written = []
    for sh in shapes:
        sid = sh["shape_id"]
        placements = [p for p in places if p["shape_id"] == sid]
        ref = placements[0]["solid_id"] if placements else None
        bores = [b for b in bores_all if b["solid_id"] == ref]
        geom = S.prism_geometry([r for r in prisms_all
                                 if r["shape_id"] == sid])

        pt, why = placement_point(sh, placements)

        # Reuse the single-shape renderer, then add the origin marker to
        # whatever 3D axes it built.
        S.draw(sh, bores, placements, None, None, geom)
        fig = plt.gcf()
        # The PART view, not merely the first 3D axes. show_shape.draw()
        # creates the context inset first, so ax3d[0] is that inset -- which
        # is where every crosshair silently went after the grid layout landed,
        # invisible against the model wireframe. draw() tags the panels it
        # reserves, so identify the part view by exclusion.
        # Every panel draw() reserves, not just two. The "as placed" pad is
        # also a 3D axes, and leaving it out of this set sent the crosshair
        # there instead of to the part view -- the same regression as when the
        # context inset was the only one excluded.
        reserved = {id(getattr(fig, "_stm_context", None)),
                    id(getattr(fig, "_stm_text", None)),
                    id(getattr(fig, "_stm_placed", None))}
        part_axes = [a for a in fig.axes
                     if hasattr(a, "get_zlim") and id(a) not in reserved]
        if part_axes:
            a = part_axes[0]
            mark_origin(a, (a.get_xlim(), a.get_ylim(), a.get_zlim()))
        # Text goes in the reserved lower-left cell, under the context view,
        # rather than floating at the figure corner: with a fixed canvas the
        # block has a known amount of room, so it can list more placements
        # without running into a panel.
        # Draw the rotated view first: it picks the representative copy, and
        # the context view then highlights that same solid in the same colour.
        chosen = add_placed(fig, sh, placements, geom, bores)
        add_context(fig, placements, all_boxes,
                    chosen["solid_id"] if chosen else None)
        ax_text = getattr(fig, "_stm_text", None)
        if pt is not None and ax_text is not None:
            ax_text.text(0.0, 1.0,
                         "origin (red crosshair) = %s\n\n%s\n\n%s"
                         % (why, placement_lines(placements, limit=9),
                            g4_excerpt(sh, placements, chosen)),
                         fontsize=6.4, va="top", ha="left", color="#333333",
                         family="monospace", transform=ax_text.transAxes)

        fname = "shape_%02d_%s_%s.png" % (int(sid), sh["type"],
                                          safe(sh.get("name")))
        path = os.path.join(args.outdir, fname)
        fig.savefig(path, dpi=args.dpi)
        plt.close(fig)
        written.append(fname)
        print("  %s" % fname)

    # The model-wide tilt figure belongs in the same directory: it is the only
    # plot that shows all 224 solids at once, and it is what says which of them
    # the de-tilt applies to.
    tilt_png = os.path.join(args.outdir, "_model_tilt_state.png")
    rc = subprocess.call([sys.executable,
                          os.path.join(HERE, "plot_tilt_state.py"),
                          "--save", tilt_png])
    if rc == 0:
        written.append(os.path.basename(tilt_png))
        print("  %s" % os.path.basename(tilt_png))
    else:
        print("  WARNING: plot_tilt_state.py failed (exit %d)" % rc)

    print("\n%d files in %s" % (len(written), args.outdir))
    return 0


if __name__ == "__main__":
    sys.exit(main())
