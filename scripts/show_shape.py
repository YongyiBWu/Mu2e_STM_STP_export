"""
Draw a shape or a placed solid from the extractor's CSVs, with dimensions
and material.

    "C:/Program Files/FreeCAD 1.1/bin/python.exe" scripts/show_shape.py --shape 1
    "C:/Program Files/FreeCAD 1.1/bin/python.exe" scripts/show_shape.py --solid 47
    "C:/Program Files/FreeCAD 1.1/bin/python.exe" scripts/show_shape.py --shape 2 --save b.png
    "C:/Program Files/FreeCAD 1.1/bin/python.exe" scripts/show_shape.py --list

    --shape N     a row of stm_shapes.csv: the block definition
    --solid N     a row of stm_placements.csv: draws its shape, plus where it sits
    --list        table of every shape, no plot
    --save FILE   write a PNG instead of opening a window
    --no-show     compute and report, draw nothing (for checking a shape quickly)

What is drawn
    BOX        the block, edges to scale.
    BOX_HOLE   the block, plus each bore as a cylinder through it, at its real
               radius and position, along its own axis.
    TUBE       a hollow cylinder, inner and outer radius, along its axis.
    PRISM      the ENVELOPE only, drawn translucent and labelled as such.

Read the PRISM case carefully
    The extractor does not model prisms: they are wedges, chamfered and skewed
    blocks that need G4Trap or G4ExtrudedSolid, and it reports only their
    bounding envelope (see the PRISM branch of classify() in
    extract_stm_geometry.py). Drawing one as a solid box would invent geometry
    that is not in the CSV, so the envelope is drawn as a translucent wireframe
    with a warning in the title. The real shape is SMALLER than what you see.

    Likewise rotY45=yes means the solid keeps off-axis faces after the de-tilt
    -- it is turned relative to the Mu2e axes. The CSV carries the flag but not
    the rotation, so the drawing shows the block in its own frame, axis-aligned,
    and says so. Do not read an orientation off these plots for such shapes.

Everything else is honest
    For BOX, BOX_HOLE and TUBE the dx/dy/dz are true world-frame extents, so
    those drawings are to scale and correctly oriented.
"""

import argparse
import collections
import csv
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
OUTDIR = os.path.join(ROOT, "output")

# CAD colour name -> something plausible to draw it in. Only for the picture;
# the material column is what carries meaning.
DRAW_COLOUR = {
    "granite gray": "#6e7b8b",
    "medium steel": "#8c9aa6",
    "medium maroon": "#b06a4e",
    "strong ice": "#9fc7d8",
    "silver gray": "#c0c4c8",
    "medium ice": "#cfe3ec",
    "pale sky": "#b9d3e8",
    "powder gray": "#9a9a9a",
    "ash gray": "#a9a9a9",
    "pale red": "#d08f8f",
}
DEFAULT_COLOUR = "#95a5a6"


def load(name, required=True):
    path = os.path.join(OUTDIR, name)
    if not os.path.exists(path):
        if not required:
            return []
        sys.exit("missing %s -- run scripts/extract_stm_geometry.py first" % path)
    with open(path) as fh:
        return list(csv.DictReader(fh))


def prism_geometry(rows):
    """Ordered cap outline, sweep axis and length, from stm_prisms.csv rows.

    Returns None when the shape has no outline, in which case only the
    bounding envelope is known and the caller must say so.
    """
    if not rows:
        return None
    rows = sorted(rows, key=lambda r: int(r["seq"]))
    r0 = rows[0]
    return {
        "outline": [(num(r, "u", 0.0), num(r, "v", 0.0)) for r in rows],
        "len": num(r0, "len", 0.0),
        "area": num(r0, "area", 0.0),
        "axis": (num(r0, "axis_x", 0.0), num(r0, "axis_y", 0.0),
                 num(r0, "axis_z", 1.0)),
        "u": (num(r0, "u_x", 1.0), num(r0, "u_y", 0.0), num(r0, "u_z", 0.0)),
        "w": (num(r0, "w_x", 0.0), num(r0, "w_y", 1.0), num(r0, "w_z", 0.0)),
    }


def num(row, key, default=None):
    """A CSV cell as float, or default when the column is blank."""
    v = (row.get(key) or "").strip()
    if not v:
        return default
    try:
        return float(v)
    except ValueError:
        return default


def mu2e_axes(x, y, z):
    """Map Mu2e (x, y, z) onto matplotlib's (horizontal, depth, vertical).

    matplotlib's 3D third argument is always the vertical axis, and no
    view_init changes that -- rotating the camera tilts the scene but leaves
    that slot upright. So +y up means passing y third and z second.

    Every function that emits 3D geometry returns points through this, so the
    whole figure shares one convention.
    """
    return x, z, y


def box_faces(cx, cy, cz, dx, dy, dz):
    """The six faces of an axis-aligned box, as lists of corner points."""
    hx, hy, hz = dx / 2.0, dy / 2.0, dz / 2.0
    x0, x1 = cx - hx, cx + hx
    y0, y1 = cy - hy, cy + hy
    z0, z1 = cz - hz, cz + hz
    faces = [
        [(x0, y0, z0), (x1, y0, z0), (x1, y1, z0), (x0, y1, z0)],   # bottom
        [(x0, y0, z1), (x1, y0, z1), (x1, y1, z1), (x0, y1, z1)],   # top
        [(x0, y0, z0), (x1, y0, z0), (x1, y0, z1), (x0, y0, z1)],   # front
        [(x0, y1, z0), (x1, y1, z0), (x1, y1, z1), (x0, y1, z1)],   # back
        [(x0, y0, z0), (x0, y1, z0), (x0, y1, z1), (x0, y0, z1)],   # left
        [(x1, y0, z0), (x1, y1, z0), (x1, y1, z1), (x1, y0, z1)],   # right
    ]
    return [[mu2e_axes(*p) for p in f] for f in faces]


def basis(axis):
    """Two unit vectors square to `axis`, for sweeping a circle around it."""
    ax = list(axis)
    n = math.sqrt(sum(v * v for v in ax)) or 1.0
    ax = [v / n for v in ax]
    # Pick whichever world axis is least parallel, so the cross product is sane.
    other = [1.0, 0.0, 0.0]
    if abs(ax[0]) > 0.9:
        other = [0.0, 1.0, 0.0]
    u = [ax[1] * other[2] - ax[2] * other[1],
         ax[2] * other[0] - ax[0] * other[2],
         ax[0] * other[1] - ax[1] * other[0]]
    un = math.sqrt(sum(v * v for v in u)) or 1.0
    u = [v / un for v in u]
    w = [ax[1] * u[2] - ax[2] * u[1],
         ax[2] * u[0] - ax[0] * u[2],
         ax[0] * u[1] - ax[1] * u[0]]
    return ax, u, w


def cylinder_faces(centre, axis, radius, length, seg=24):
    """Quad strips forming the wall of a cylinder."""
    ax, u, w = basis(axis)
    cx, cy, cz = centre
    h = length / 2.0
    ring = []
    for i in range(seg):
        t = 2.0 * math.pi * i / seg
        off = [radius * (math.cos(t) * u[k] + math.sin(t) * w[k]) for k in range(3)]
        ring.append(off)
    faces = []
    for i in range(seg):
        a, b = ring[i], ring[(i + 1) % seg]
        p0 = (cx + a[0] - ax[0] * h, cy + a[1] - ax[1] * h, cz + a[2] - ax[2] * h)
        p1 = (cx + b[0] - ax[0] * h, cy + b[1] - ax[1] * h, cz + b[2] - ax[2] * h)
        p2 = (cx + b[0] + ax[0] * h, cy + b[1] + ax[1] * h, cz + b[2] + ax[2] * h)
        p3 = (cx + a[0] + ax[0] * h, cy + a[1] + ax[1] * h, cz + a[2] + ax[2] * h)
        faces.append([mu2e_axes(*p) for p in (p0, p1, p2, p3)])
    return faces


def describe(shape, bores, placements, geom=None):
    """The text half of the answer: what this shape is, in words."""
    sid = shape["shape_id"]
    kind = shape["type"]
    dx = num(shape, "dx", 0.0)
    dy = num(shape, "dy", 0.0)
    dz = num(shape, "dz", 0.0)
    lines = []
    lines.append("shape %s  %s  %s" % (sid, kind, shape.get("name") or "(unnamed)"))
    lines.append("  placed %s time(s)" % shape.get("count"))
    # A prism has no dx/dy/dz: it is a swept outline, not a box, and the
    # extractor blanks those columns rather than publish an envelope that
    # could be read as a buildable size.
    if None not in (num(shape, "dx"), num(shape, "dy"), num(shape, "dz")):
        lines.append("  size      %.3f x %.3f x %.3f mm" % (dx, dy, dz))
        ix = num(shape, "in_x")
        iy = num(shape, "in_y")
        iz = num(shape, "in_z")
        if None not in (ix, iy, iz):
            lines.append("            %.3f x %.3f x %.3f in" % (ix, iy, iz))
    lines.append("  material  %s" % (shape.get("material") or "?"))
    hint = shape.get("material_hint")
    if hint:
        lines.append("  hint      %s   (a guess, not authoritative)" % hint)
    src = shape.get("material_src")
    if src:
        lines.append("  from      %s" % src)

    if kind == "TUBE":
        lines.append("  bore      rmin %.3f  rmax %.3f mm"
                     % (num(shape, "rmin", 0.0), num(shape, "rmax", 0.0)))
    if bores:
        lines.append("  %d bore(s):" % len(bores))
        for b in bores:
            lines.append("      r=%.3f at (%.2f, %.2f, %.2f) along (%.3f, %.3f, %.3f)"
                         % (num(b, "r", 0.0), num(b, "x", 0.0), num(b, "y", 0.0),
                            num(b, "z", 0.0), num(b, "ax", 0.0),
                            num(b, "ay", 0.0), num(b, "az", 0.0)))

    if geom is not None:
        lines.append("  extrusion %d-sided cap, area %.1f mm^2, swept %.3f mm"
                     % (len(geom["outline"]), geom["area"], geom["len"]))
        lines.append("            axis (%.4f, %.4f, %.4f)" % geom["axis"])
        lines.append("            -> G4ExtrudedSolid; outline in stm_prisms.csv")

    warn = []
    if kind == "PRISM" and geom is None:
        warn.append("PRISM: the CSV holds only the bounding ENVELOPE. The real "
                    "solid is a wedge or chamfered block and is SMALLER than "
                    "what is drawn.")
    if kind == "OTHER":
        warn.append("OTHER: unclassified; only the envelope is known.")
    if (shape.get("rotY45") or "").strip():
        if geom is not None:
            # For an extrusion the orientation IS recoverable: stm_prisms.csv
            # carries the sweep axis and the cap plane's basis in Mu2e
            # coordinates. The flag only means the outline is cut at 45 deg
            # within its own plane, which the drawing shows correctly.
            warn.append("rotY45: the cap outline has 45 deg features. The "
                        "orientation is NOT lost -- axis_*/u_*/w_* in "
                        "stm_prisms.csv place this solid in Mu2e coordinates.")
        else:
            warn.append("rotY45: this solid is turned relative to the Mu2e "
                        "axes. The CSV carries the flag but not the rotation, "
                        "so it is drawn in its own frame. Do not read "
                        "orientation off this.")
    mm = [num(p, "material_match") for p in placements]
    mm = [v for v in mm if v is not None]
    if mm and min(mm) < 0.9:
        warn.append("material_match as low as %.2f: the material came from a "
                    "loose geometric match and may be wrong." % min(mm))
    return lines, warn


def prism_faces(geom):
    """The swept solid: two caps plus one quad per outline edge.

    The outline is 2D in the cap's own plane, so each point is lifted into 3D
    and offset +/- half the length along the sweep axis.

    Drawn EXACTLY as stored -- the outline is not re-centred. The extractor
    anchors it on vertex 0, so local (0,0) is a corner of the polygon, and
    re-centring here moved the drawing off that corner: the origin crosshair
    landed 349 mm from the nearest vertex on shape 3 and 137 mm on shape 23.
    A box is centred on its origin and a prism is not; that difference is real
    and the picture should show it.
    """
    half = geom["len"] / 2.0
    pts = geom["outline"]

    def lift(p, s):
        # LOCAL frame: the cap lies in x,y and the sweep runs along z, exactly
        # as G4ExtrudedSolid receives it. This deliberately does NOT apply the
        # stored u/w/axis basis.
        #
        # Lifting through that basis here produced world-placed geometry, which
        # broke the convention twice: the part view is meant to show the shape
        # in its own frame before any rotation, and nothing else then drew the
        # local frame at all. The orientation belongs to the placement, and the
        # "as placed" pad applies it.
        return (p[0], p[1], s * half)

    lo = [lift(p, -1.0) for p in pts]
    hi = [lift(p, +1.0) for p in pts]
    faces = [lo, hi]
    for i in range(len(pts)):
        j = (i + 1) % len(pts)
        faces.append([lo[i], lo[j], hi[j], hi[i]])
    return [[mu2e_axes(*p) for p in f] for f in faces]


def draw(shape, bores, placements, save, solid_row=None, geom=None):
    import matplotlib
    if save:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection

    kind = shape["type"]
    dx = num(shape, "dx", 0.0)
    dy = num(shape, "dy", 0.0)
    dz = num(shape, "dz", 0.0)
    colour = DRAW_COLOUR.get((shape.get("material") or "").strip(), DEFAULT_COLOUR)
    # A prism with an outline is drawn exactly, not as an envelope: the CSV
    # carries the real swept polygon, so claiming "envelope only" would now be
    # a lie about data we actually have.
    envelope = kind in ("PRISM", "OTHER", "PRISM_HOLE") and geom is None

    # Extreme aspect ratios are the norm here: shape 2 is 914 x 20 x 317, a
    # 46:1 slab. A lone 3D view of that is nearly unreadable -- the short axis
    # collapses and its tick labels overprint -- and a 2.5mm bore in a 914mm
    # face is a third of a percent of the frame. So pair the 3D view with a
    # flat cross-section through the largest face, where bores are locatable.
    # A bored shape gets two sections, not one: across the bore axis, where a
    # hole is a circle and its radius and position read directly, and ALONG it,
    # where the hole is a rectangle showing how deep it goes and whether it
    # passes through. Neither view alone answers both questions -- the
    # perpendicular one hides depth, the parallel one hides radius.
    flat = kind != "TUBE" and bool(bores)

    # ONE canvas size for every shape, so a directory of these can be flipped
    # through without each page jumping to a different scale. The layout is a
    # fixed 2-column grid:
    #
    #   left  column : context (where this part sits in the model), then the
    #                  text block underneath it
    #   right column : the part itself, then the section views if it has bores
    #
    # Shapes without bores leave the lower-right cells empty rather than
    # resizing the figure -- an empty cell costs nothing and a moving canvas
    # costs comparability.
    fig = plt.figure(figsize=(16.0, 9.0))
    # The context view needs most of the left column's height -- it is a 3D
    # plot of the whole model and shrinks badly -- while the text block under
    # it is a short monospace list. Weighting the rows 1.45:1 gives the
    # context room without leaving the text stranded at the bottom.
    # Six rows so the left column can hold three stacked panels while the
    # right-hand panels still split cleanly in half (3 rows each).
    gs = fig.add_gridspec(6, 3, width_ratios=[1.0, 1.25, 1.25],
                          left=0.02, right=0.98, top=0.92, bottom=0.05,
                          wspace=0.24, hspace=0.22)
    # Left column, top to bottom: where the part sits in the model, the same
    # part AFTER its rotation, then the text. The middle pad is what makes the
    # rotation legible -- the main panel necessarily draws the shape in its
    # CANONICAL frame (that is what the dimensions describe), so without a
    # rotated view nothing on the page connects those dimensions to the
    # orientation column.
    ax_context = fig.add_subplot(gs[0:2, 0], projection="3d")
    ax_placed = fig.add_subplot(gs[2:4, 0], projection="3d")
    ax_text = fig.add_subplot(gs[4:6, 0])
    ax_text.axis("off")
    fig._stm_context = ax_context
    fig._stm_placed = ax_placed
    fig._stm_text = ax_text

    # The part view spans both rows of column 2 when there are no sections to
    # show, so a plain box still gets a large readable panel.
    if flat:
        ax = fig.add_subplot(gs[:, 1], projection="3d")
    else:
        ax = fig.add_subplot(gs[:, 1:], projection="3d")

    if kind == "TUBE":
        rmax = num(shape, "rmax", max(dx, dy) / 2.0)
        rmin = num(shape, "rmin", 0.0)
        # The CSV gives the tube's length in dz and its axis is not stored per
        # shape, so draw along z: for both tubes here that is the real axis.
        for r, alpha in ((rmax, 0.55), (rmin, 0.9)):
            if r and r > 0:
                fs = cylinder_faces((0.0, 0.0, 0.0), (0, 0, 1), r, dz)
                pc = Poly3DCollection(fs, alpha=alpha, facecolor=colour,
                                      edgecolor="none")
                ax.add_collection3d(pc)
        lim = [(-rmax, rmax), (-rmax, rmax), (-dz / 2.0, dz / 2.0)]
        aspect = (2 * rmax, 2 * rmax, dz)
    elif geom is not None:
        fs = prism_faces(geom)
        ax.add_collection3d(Poly3DCollection(
            fs, alpha=0.45, facecolor=colour, edgecolor="#2c3e50",
            linewidth=1.0))
        # prism_faces() already emits DISPLAY order (Mu2e x, z, y) via
        # mu2e_axes. The shared limit code below applies the swap itself --
        # set_ylim(lim[2]), set_zlim(lim[1]) -- so it expects lim in MU2E
        # order, the way the box branch supplies it from dx/dy/dz.
        #
        # Deriving lim from the swapped points and letting it be swapped again
        # put the 150.4 mm sweep on the axis labelled "y [+y up]" and made the
        # local frame look like the placed one. Convert back to Mu2e order here
        # so both branches speak the same convention.
        pts = [p for f in fs for p in f]
        disp = [(min(p[k] for p in pts), max(p[k] for p in pts))
                for k in range(3)]
        lim = [disp[0], disp[2], disp[1]]          # (x, z, y) -> (x, y, z)
        aspect = tuple(max(hi - lo, 1e-6) for lo, hi in lim)
    else:
        fs = box_faces(0.0, 0.0, 0.0, dx, dy, dz)
        # A literal "none" facecolor gives mplot3d an empty colour array, and
        # its do_3d_projection then unpacks five values from nothing and
        # raises. A fully transparent RGBA is the same picture and survives.
        face = (0.0, 0.0, 0.0, 0.0) if envelope else colour
        pc = Poly3DCollection(
            fs,
            alpha=0.18 if envelope else 0.30,
            facecolor=face,
            edgecolor="#c0392b" if envelope else "#2c3e50",
            linewidth=1.6 if envelope else 1.0,
            linestyle="--" if envelope else "-")
        ax.add_collection3d(pc)

        # Bores are given in Mu2e coordinates; the block is drawn centred on the
        # origin, so shift each bore by the solid's own centre to line them up.
        if bores:
            for b in bores:
                r = num(b, "r", 0.0)
                # dx/dy/dz are already relative to the block centre, which is
                # what the block is drawn around. The absolute x/y/z in the CSV
                # are for cross-checking against CAD, not for drawing.
                bx = num(b, "dx", 0.0)
                by = num(b, "dy", 0.0)
                bz = num(b, "dz", 0.0)
                axis = (num(b, "ax", 0.0), num(b, "ay", 0.0), num(b, "az", 1.0))
                # Run the bore right through the block on its own axis. Drawn
                # at TRUE radius -- never enlarged to be visible, since this
                # tool exists to report dimensions faithfully. The 2D panel is
                # what makes a small bore findable.
                # Draw the bore its OWN depth, not the block's full extent.
                # Using the extent makes every bore look like a through-hole:
                # shape 8's r=12.7 counterbore is 50.602 mm deep in a 101.6 mm
                # block, so it was drawn at twice its length and appeared to
                # run right through. The 2D parallel section already used
                # depth, which is why the two views disagreed.
                depth = num(b, "depth", 0.0)
                if depth <= 0:
                    # No depth recorded: fall back to the block extent along
                    # the bore axis, which is the old behaviour and right for
                    # a genuine through-hole.
                    depth = (abs(axis[0]) * dx + abs(axis[1]) * dy
                             + abs(axis[2]) * dz)
                fs = cylinder_faces((bx, by, bz), axis, r, depth)
                ax.add_collection3d(Poly3DCollection(
                    fs, alpha=1.0, facecolor="#c0392b", edgecolor="#7b241c",
                    linewidth=0.3))
        lim = [(-dx / 2.0, dx / 2.0), (-dy / 2.0, dy / 2.0), (-dz / 2.0, dz / 2.0)]
        aspect = (dx, dy, dz)

    # Pad slightly so edges are not flush against the axes box.
    lim = [(lo - (hi - lo) * 0.08, hi + (hi - lo) * 0.08) for lo, hi in lim]
    # lim and aspect are in Mu2e (x, y, z); the axes hold (x, z, y).
    ax.set_xlim(*lim[0])
    ax.set_ylim(*lim[2])
    ax.set_zlim(*lim[1])
    try:
        # Aspect proportional to the real extents, so a 914 x 20 x 317 slab is
        # drawn as a slab. A cube box aspect (the obvious choice) would show
        # every shape as a cube-ish blob and hide exactly the proportions the
        # plot exists to convey.
        ax.set_box_aspect((aspect[0], aspect[2], aspect[1]))
    except (AttributeError, ValueError):
        pass
    # Mu2e display convention: +y up, -x to the right, +z into the page. The
    # data swap lives in mu2e_axes(); here the horizontal axis is inverted so
    # -x runs right, and each axis is named for the Mu2e quantity on it rather
    # than for matplotlib's letter for that slot.
    # Raised viewpoint: looking down from well above +y reads the model far
    # better than a near-horizontal view, where the shield house's own walls
    # hide everything behind them.
    ax.view_init(elev=32, azim=-68)
    ax.set_xlim(ax.get_xlim()[::-1])
    ax.set_xlabel("x (mm)  [-x right]")
    ax.set_ylabel("z (mm)  [+z into page]")
    ax.set_zlabel("y (mm)  [+y up]")
    # On a very flat solid the short axis's ticks collide into an unreadable
    # smear, so thin them out. The dimensions are in the title anyway.
    for axis, extent in zip((ax.xaxis, ax.yaxis, ax.zaxis), aspect):
        if extent < max(aspect) * 0.25:
            axis.set_major_locator(plt.MaxNLocator(3))

    title = "shape %s  %s  %s" % (shape["shape_id"], kind,
                                  shape.get("name") or "")
    sub = "%.1f x %.1f x %.1f mm   %s" % (dx, dy, dz,
                                          shape.get("material") or "?")
    if kind == "TUBE":
        sub = "rmin %.2f  rmax %.2f  length %.1f mm   %s" % (
            num(shape, "rmin", 0.0), num(shape, "rmax", 0.0), dz,
            shape.get("material") or "?")
    if solid_row is not None:
        sub += "\nsolid %s at (%.1f, %.1f, %.1f) mm" % (
            solid_row["solid_id"], num(solid_row, "x", 0.0),
            num(solid_row, "y", 0.0), num(solid_row, "z", 0.0))
    if envelope:
        sub += "\nENVELOPE ONLY -- real solid is smaller"
    elif geom is not None:
        sub = ("%d-sided cap, area %.0f mm^2, swept %.1f mm   %s"
               % (len(geom["outline"]), geom["area"], geom["len"],
                  shape.get("material") or "?"))
        sub += "\nexact extrusion from stm_prisms.csv"
    if (shape.get("rotY45") or "").strip() and geom is None:
        sub += "\nrotY45: drawn in its own frame, orientation not shown"
    ax.set_title(title + "\n" + sub, fontsize=10)

    if flat:
        # Look down the block's shortest axis: that is the face the bores go
        # through, and the one with room to show them.
        ext = [dx, dy, dz]
        # A prism carries no dx/dy/dz -- stm_shapes.csv leaves them blank for
        # every swept outline, which is the right call for a shape that is a
        # polygon and a length rather than a box. Until PRISM_HOLE there was no
        # shape that was BOTH a prism and bored, so nothing reached this
        # section code with empty dims; shape 24 did, and ext = [0,0,0] divided
        # by zero at the panel-alignment step.
        #
        # The outline knows the answer: its frame is the canonical one (cap on
        # x,y, sweep on z), so the spans in u and v plus the sweep length ARE
        # the extents the sections need.
        if geom is not None and not any(ext):
            us = [p[0] for p in geom["outline"]]
            vs = [p[1] for p in geom["outline"]]
            ext = [max(us) - min(us), max(vs) - min(vs), geom["len"]]

        # Lower corner of the block in each canonical axis.
        #
        # Everything below used to assume the block straddles the origin, and
        # drew it as a rectangle from -ext/2. That holds for a box, whose bore
        # offsets are measured from its centre, but NOT for a prism: the
        # outline is anchored on a VERTEX, so it spans min(u)..max(u) rather
        # than +/-ext/2, and its bore offsets are measured from that same
        # anchor. Drawing shape 29's outline (u 0..203.2, v 0..50.8) centred
        # put its bore at (105.4, 25.4) on axes running +/-101.6 and +/-25.4 --
        # out at the corner, half of it outside the block, which read both as
        # a misplaced hole and as one that failed to cut through.
        #
        # Only the cap plane is anchored; the sweep still straddles the origin
        # (-len/2..+len/2), which is what G4ExtrudedSolid builds and what the
        # bore's own dz is measured against.
        if geom is not None:
            us = [p[0] for p in geom["outline"]]
            vs = [p[1] for p in geom["outline"]]
            lo = [min(us), min(vs), -ext[2] / 2.0]
        else:
            lo = [-ext[k] / 2.0 for k in range(3)]
        # Section PERPENDICULAR TO THE BORE AXIS, not across the thinnest
        # extent. Choosing the thinnest axis sections the plane that CONTAINS
        # the bore axes, where a hole is a rectangle seen edge-on -- drawing it
        # as a circle there invents geometry and makes coaxial bores at
        # different depths look like separate holes side by side. Shape 8 is
        # the case that exposed it: two z-running bores counterbored r=12.7 to
        # r=6.985 appeared as four parallel holes.
        axis_votes = collections.Counter()
        for b in bores:
            a = (abs(num(b, "ax", 0.0)), abs(num(b, "ay", 0.0)),
                 abs(num(b, "az", 0.0)))
            axis_votes[max(range(3), key=lambda i: a[i])] += 1
        if axis_votes:
            thin = axis_votes.most_common(1)[0][0]
        else:
            thin = min(range(3), key=lambda i: ext[i])
        horiz, vert = [i for i in range(3) if i != thin]
        names = ["x", "y", "z"]
        # (1, 3, 2), not (1, 2, 2): the figure is a three-column grid, so
        # asking for "the second of two" claims the whole right half and the
        # third panel then lands on top of it.
        # Column 3, top: the section across the bore axis.
        # Rows 0:3 of the six-row grid -- the top half. Indexing gs[0, 2] here
        # would claim a single row, a sixth of the height, now that the left
        # column needs six rows to stack three panels.
        # A prism's cap is a POLYGON, not a rectangle, and exactly one of the
        # two section panels looks down on that cap. Which one is not fixed:
        # the panels are chosen from the bore axis, so for shapes 24 and 29 --
        # bored across the cap rather than along the sweep -- the cap plane
        # landed on panel 3, not panel 2. Hardcoding either would draw the
        # outline on the wrong panel for a prism bored the other way.
        #
        # The cap plane is canonical (u,v) = axes 0 and 1, so a panel shows it
        # iff its two axes are exactly {0, 1}.
        def block_patch(h, v, **kw):
            """The block's outline in the (h, v) plane."""
            if geom is not None and {h, v} == {0, 1}:
                pts = [(p[0], p[1]) if h == 0 else (p[1], p[0])
                       for p in geom["outline"]]
                return plt.Polygon(pts, closed=True, **kw)
            return plt.Rectangle((lo[h], lo[v]), ext[h], ext[v], **kw)

        ax2 = fig.add_subplot(gs[0:3, 2])
        ax2.add_patch(block_patch(
            horiz, vert,
            facecolor=colour if not envelope else "none", alpha=0.30,
            edgecolor="#c0392b" if envelope else "#2c3e50",
            linestyle="--" if envelope else "-", linewidth=1.4))
        # dx/dy/dz are ALREADY relative to the block centre, which the section
        # is drawn around. Subtracting the placement again (a leftover from the
        # absolute-coordinate schema) shifted every bore by the solid's own
        # position.
        #
        # Draw the widest bore first so a counterbore nests visibly inside its
        # pilot hole rather than hiding it.
        order = sorted(bores, key=lambda b: -num(b, "r", 0.0))
        seen = []
        for b in order:
            p = (num(b, "dx", 0.0), num(b, "dy", 0.0), num(b, "dz", 0.0))
            r = num(b, "r", 0.0)
            c = (p[horiz], p[vert])
            ax2.add_patch(plt.Circle(c, r, facecolor="#c0392b",
                                     edgecolor="#7b241c", alpha=0.75,
                                     linewidth=0.8, zorder=3))
            # Two bores projecting onto the same point are coaxial: label them
            # once, stacked, instead of printing overlapping text.
            dup = [q for q in seen if abs(q[0] - c[0]) < 1e-3
                   and abs(q[1] - c[1]) < 1e-3]
            ax2.annotate("r%.2f" % r, c, textcoords="offset points",
                         xytext=(6, 6 + 11 * len(dup)), fontsize=7.5,
                         color="#7b241c", zorder=4)
            seen.append(c)
        # set_aspect("equal") overrides set_ylim, so set the aspect FIRST and
        # then pad the short axis by hand. Doing it the other way round leaves
        # matplotlib's expanded range and its ticks printed over ours.
        # Set the limits first, then lock the aspect with adjustable="box":
        # that keeps the data range we asked for and reshapes the axes instead.
        # With the default adjustable="datalim" matplotlib widens the short
        # axis to match the long one, which on a 406 x 102 block left the
        # section floating in whitespace at +/-200.
        ax2.set_xlim(lo[horiz] - ext[horiz] * 0.04,
                     lo[horiz] + ext[horiz] * 1.04)
        ax2.set_ylim(lo[vert] - ext[vert] * 0.04,
                     lo[vert] + ext[vert] * 1.04)
        ax2.set_aspect("equal", adjustable="box")
        ax2.set_xlabel("%s (mm)" % names[horiz])
        ax2.set_ylabel("%s (mm)" % names[vert])

        # --- third panel: the section ALONG the bore axis ---------------
        # Same block, cut on the plane the bores run in, so each hole shows as
        # a rectangle: its length is how deep the bore goes and whether it
        # breaks through the far face. The perpendicular panel cannot show
        # that, and a counterbore looks identical to a through-hole there.
        # Column 3, bottom: the section along the bore axis, directly under the
        # perpendicular one so the two read as a pair.
        ax3 = fig.add_subplot(gs[3:6, 2])
        long_h = horiz            # keep the same horizontal axis as panel 2
        long_v = thin             # vertical axis is now the bore direction
        ax3.add_patch(block_patch(
            long_h, long_v,
            facecolor=colour if not envelope else "none", alpha=0.30,
            edgecolor="#c0392b" if envelope else "#2c3e50",
            linestyle="--" if envelope else "-", linewidth=1.4))
        for b in sorted(bores, key=lambda b: -num(b, "r", 0.0)):
            p = (num(b, "dx", 0.0), num(b, "dy", 0.0), num(b, "dz", 0.0))
            r = num(b, "r", 0.0)
            depth = num(b, "depth", 0.0)
            if depth <= 0:
                depth = ext[long_v]
            # The bore is centred on its own midpoint along the axis, which is
            # what the depth column measures.
            ax3.add_patch(plt.Rectangle(
                (p[long_h] - r, p[long_v] - depth / 2.0), 2 * r, depth,
                facecolor="#c0392b", edgecolor="#7b241c", alpha=0.75,
                linewidth=0.8, zorder=3))
        ax3.set_xlim(lo[long_h] - ext[long_h] * 0.04,
                     lo[long_h] + ext[long_h] * 1.04)
        ax3.set_ylim(lo[long_v] - ext[long_v] * 0.04,
                     lo[long_v] + ext[long_v] * 1.04)
        # Deliberately NOT adjustable="box" here. That mode lets matplotlib
        # resize the axes box to suit the data, and it recomputes at draw time,
        # which silently undoes the explicit alignment applied below. The box
        # is instead sized by hand, and "equal" scale is achieved by giving it
        # a height proportional to its own y range at panel 2's x scale.
        ax3.set_aspect("auto")
        ax3.set_xlabel("%s (mm)" % names[long_h])
        ax3.set_ylabel("%s (mm)" % names[long_v])
        ax3.set_title("section through %s%s, along the bore axis (%s)\n"
                      "bore depth and break-through, true scale"
                      % (names[long_h], names[long_v], names[thin]),
                      fontsize=9)
        ax3.grid(alpha=0.25, linewidth=0.5)
        ax2.set_title("section through %s%s, looking along %s\n"
                      "%d bore(s), true scale"
                      % (names[horiz], names[vert], names[thin], len(bores)),
                      fontsize=9)
        ax2.grid(alpha=0.25, linewidth=0.5)

        # Align the two sections: same width, same left edge.
        #
        # Both already span the same x range -- long_h IS horiz -- but
        # set_aspect("equal", adjustable="box") shrinks each axes box to suit
        # its OWN data aspect, and their vertical extents differ (81.28 vs
        # 101.6 on shape 8). That left the panels different widths and
        # horizontally offset, so the same x coordinate sat at two different
        # places on the page and the views could not be read against each
        # other.
        #
        # Fix: let matplotlib settle panel 2 (which keeps adjustable="box", so
        # its own aspect is honest), then give panel 3 that exact x-position
        # and width, and compute its HEIGHT so that millimetres-per-inch match
        # panel 2 horizontally and vertically. That reproduces equal aspect by
        # construction instead of asking matplotlib for it, so nothing is
        # recomputed behind us at save time.
        fig.canvas.draw()
        p2 = ax2.get_position()
        p3 = ax3.get_position()
        fw, fh = fig.get_size_inches()

        # Size BOTH panels from one shared width, chosen so that each keeps
        # true scale (equal mm per inch in x and y) AND fits its own cell.
        #
        # Letting matplotlib pick per-panel widths is what broke alignment in
        # the first place; picking a width for panel 2 alone and forcing it on
        # panel 3 fixed alignment but left panel 3 slightly squashed when its
        # true-scale height exceeded the cell -- 2.7% on shape 8. Deriving one
        # width from the tighter of the two constraints makes both exact.
        x_mm = ext[long_h] * 1.08               # the +/-0.54 range, doubled
        v2_mm = ext[vert] * 1.08                # panel 2's vertical span
        v3_mm = ext[long_v] * 1.08              # panel 3's vertical span

        # Widest this panel may be before its true-scale height overflows the
        # cell, per panel; take the smaller so neither overflows.
        cap2 = (p2.height * fh) * (x_mm / v2_mm) / fw
        cap3 = (p3.height * fh) * (x_mm / v3_mm) / fw
        width = min(p2.width, cap2, cap3)

        # True-scale heights at that shared width.
        in_per_mm = (width * fw) / x_mm
        h2 = (v2_mm * in_per_mm) / fh
        h3 = (v3_mm * in_per_mm) / fh

        # Keep each panel anchored to the top of its own cell so the pair stays
        # visually joined, and keep panel 2's left edge as the common one.
        ax2.set_aspect("auto")                  # box is ours to set now
        ax2.set_position([p2.x0, p2.y0 + p2.height - h2, width, h2])
        ax3.set_position([p2.x0, p3.y0 + p3.height - h3, width, h3])

    # No tight_layout: the GridSpec above already fixes every panel's position,
    # and letting matplotlib re-flow it would make the canvas shape-dependent
    # again, which is the thing this layout exists to prevent.
    if save:
        fig.savefig(save, dpi=130)
        print("wrote %s" % save)
    else:
        plt.show()


def list_shapes(shapes):
    print("%-4s %-9s %-5s %-26s %-14s %s"
          % ("id", "type", "n", "name", "material", "size (mm)"))
    for s in shapes:
        # Prisms report their extrusion instead of a box size, since the
        # envelope columns are blank for them by design.
        if (s.get("dx") or "").strip():
            size = "%.1f x %.1f x %.1f" % (num(s, "dx", 0.0), num(s, "dy", 0.0),
                                           num(s, "dz", 0.0))
        elif (s.get("cap_area") or "").strip():
            size = "cap %.0f mm^2 x %s swept %.1f" % (
                num(s, "cap_area", 0.0), s.get("cap_verts") or "?",
                num(s, "sweep_len", 0.0))
        else:
            size = "(no dimensions)"
        print("%-4s %-9s %-5s %-26s %-14s %s%s"
              % (s["shape_id"], s["type"], s["count"],
                 (s.get("name") or "")[:26], (s.get("material") or "?")[:14],
                 size,
                 "   rotY45" if (s.get("rotY45") or "").strip() else ""))


def main():
    ap = argparse.ArgumentParser(
        description=__doc__.strip(),
        formatter_class=argparse.RawDescriptionHelpFormatter)
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--shape", type=int, help="shape_id from stm_shapes.csv")
    g.add_argument("--solid", type=int, help="solid_id from stm_placements.csv")
    g.add_argument("--list", action="store_true", help="list all shapes")
    ap.add_argument("--save", metavar="FILE", help="write a PNG instead of a window")
    ap.add_argument("--no-show", action="store_true", help="text only, no plot")
    args = ap.parse_args()

    shapes = load("stm_shapes.csv")
    if args.list:
        list_shapes(shapes)
        return 0

    places = load("stm_placements.csv")
    bores_all = load("stm_bores.csv")
    # Optional: only present once the extractor has run with prism export.
    prisms_all = load("stm_prisms.csv", required=False)

    solid_row = None
    if args.solid is not None:
        match = [p for p in places if p["solid_id"] == str(args.solid)]
        if not match:
            sys.exit("no solid_id %d in stm_placements.csv (0-%d)"
                     % (args.solid, len(places) - 1))
        solid_row = match[0]
        sid = solid_row["shape_id"]
    else:
        sid = str(args.shape)

    shape = next((s for s in shapes if s["shape_id"] == sid), None)
    if shape is None:
        sys.exit("no shape_id %s in stm_shapes.csv" % sid)

    placements = [p for p in places if p["shape_id"] == sid]
    # Bores are recorded per solid. Take them from the solid being shown, or
    # from the first placement of the shape, so the count matches nbores.
    ref_solid = solid_row["solid_id"] if solid_row else (
        placements[0]["solid_id"] if placements else None)
    bores = [b for b in bores_all if b["solid_id"] == ref_solid]

    geom = prism_geometry([r for r in prisms_all if r["shape_id"] == sid])

    lines, warn = describe(shape, bores, placements, geom)
    print("\n".join(lines))
    if solid_row is not None:
        print("  this solid: %s  at (%.2f, %.2f, %.2f) mm"
              % (solid_row.get("part") or "?", num(solid_row, "x", 0.0),
                 num(solid_row, "y", 0.0), num(solid_row, "z", 0.0)))
        mm = num(solid_row, "material_match")
        if mm is not None:
            print("  match score %.3f" % mm)
    else:
        print("  placements: %s" % ", ".join(p["solid_id"] for p in placements[:12])
              + (" ..." if len(placements) > 12 else ""))
    for w in warn:
        print("\n  ! %s" % w)

    if args.no_show:
        return 0
    try:
        draw(shape, bores, placements, args.save, solid_row, geom)
    except ImportError as exc:
        print("\n(no plot: %s)" % exc)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
