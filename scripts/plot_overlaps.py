"""Plot every 'skin' and 'real' overlap as its own figure, into output/overlap/.

    "C:/Program Files/FreeCAD 1.1/bin/python.exe" scripts/plot_overlaps.py

Why this exists
    check_overlaps.py says WHICH pairs interpenetrate and by how much. It
    cannot say whether a given pair is a designed fit or a mistake, and the
    numbers alone do not settle it: an overlap of 2007 mm^3 sounds alarming and
    is in fact a 0.778 mm eccentricity, while 9154 mm^3 is two bricks whose
    faces sit half a millimetre inside each other. Both need a picture.

    'coplanar' pairs are skipped. There are 141 of them, they are stacked
    bricks sharing a face, and nothing is learned by drawing each one.

Two families, two layouts
    The skin/real pairs split cleanly, and one layout cannot serve both:

    PLANAR   the intersection is a sheet, thin in one axis (0.51-0.86 mm here,
             always y). Two nominally-coincident faces that actually
             interpenetrate. The useful view is along the thin axis, plus a
             zoom that makes the depth legible.
    TUBULAR  a tube in a bore of the SAME radius, axes offset by a fraction of
             a millimetre. True scale shows nothing -- the two circles look
             concentric -- so the figure has to zoom onto the two axes.
    OTHER    anything matching neither; gets a generic 3D + section view rather
             than a layout that assumes structure it does not have.

A note on repetition
    The same physical fit appears in several pairs: tube 157 is reported
    against 200, 43, 162, 175, 67, 69, 220, 180 and 156, because each wall it
    passes through is a separate pair. Nine figures, one tube. That is inherent
    to a pairwise report; the per-figure text names the partner so it is clear
    which fit is being shown.
"""

import json
import math
import os
import re
import sys

import FreeCAD  # noqa: F401  (must precede Part)
import Part

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
OUTDIR = os.path.join(ROOT, "output")
PLOTDIR = os.path.join(OUTDIR, "overlap")
STP = os.path.join(ROOT, "STM_STP_files",
                   "F10269585_G4_Shield_House_3.stp")

sys.path.insert(0, HERE)
import extract_stm_geometry as E          # noqa: E402  (transform, TILT_DEG)
import show_shape as S                    # noqa: E402  (mu2e_axes, load)


# --- coordinates -------------------------------------------------------
#
# The booleans run on the STEP solids, which are in CAD coordinates. Every
# NUMBER that reaches a figure is in Mu2e coordinates instead, because these
# plots are read beside the ones in plots/ and a reader comparing them must
# not have to know that two different frames are in play. The gap is not
# subtle: solid 165's centre is x -404.254 in CAD and +49.876 in Mu2e.
#
# transform() de-tilts about z, flips 180 deg about y and shifts onto the
# anchor; transform_dir() does the rotation only, for directions.

def mu2e_pt(p):
    """A CAD point (anything with .x/.y/.z, or a 3-tuple) -> Mu2e tuple."""
    if hasattr(p, "x"):
        return E.transform(p.x, p.y, p.z)
    return E.transform(p[0], p[1], p[2])


def mu2e_dir(d):
    if hasattr(d, "x"):
        return E.transform_dir(d.x, d.y, d.z)
    return E.transform_dir(d[0], d[1], d[2])


def mu2e_bbox(shape):
    """(lo, hi) of a CAD shape's bounding box, in Mu2e coordinates.

    The transform includes a 180 deg flip, so a corner that was the minimum in
    CAD becomes the maximum in Mu2e. Transform both corners and re-sort rather
    than assuming the order survives.
    """
    b = shape.BoundBox
    pts = [E.transform(x, y, z)
           for x in (b.XMin, b.XMax)
           for y in (b.YMin, b.YMax)
           for z in (b.ZMin, b.ZMax)]
    lo = [min(p[k] for p in pts) for k in range(3)]
    hi = [max(p[k] for p in pts) for k in range(3)]
    return lo, hi


def tidy_3d(ax, pad=8.0):
    """Shrink the axis labels and tick text of a 3D panel.

    Four columns leave each 3D panel narrow enough that matplotlib's default
    label sizes run into the neighbouring panel's title. The labels still say
    which Mu2e axis is which, just smaller.
    """
    for a in (ax.xaxis, ax.yaxis, ax.zaxis):
        a.set_tick_params(labelsize=6.5, pad=1)
        a.labelpad = pad
    ax.xaxis.label.set_size(7.5)
    ax.yaxis.label.set_size(7.5)
    ax.zaxis.label.set_size(7.5)


def mu2e_3d(ax, lim):
    """Apply the Mu2e display convention to a 3D axes.

    Identical to show_shape.draw() and plot_all_components.add_context(), so
    an overlap figure and a model plot can be read against each other:
    matplotlib's third slot is always vertical, so the data goes in as
    (x, z, y) via mu2e_axes and the limits follow that order; the horizontal
    axis is then inverted so -x runs right; the viewpoint is raised.

    `lim` is a dict of Mu2e axis -> (lo, hi).
    """
    ax.set_xlim(*lim["x"])
    ax.set_ylim(*lim["z"])
    ax.set_zlim(*lim["y"])
    try:
        ax.set_box_aspect([max(lim["x"][1] - lim["x"][0], 1e-6),
                           max(lim["z"][1] - lim["z"][0], 1e-6),
                           max(lim["y"][1] - lim["y"][0], 1e-6)])
    except (AttributeError, ValueError):
        pass
    ax.view_init(elev=32, azim=-68)
    ax.set_xlim(ax.get_xlim()[::-1])
    ax.set_xlabel("x (mm)  [-x right]")
    ax.set_ylabel("z (mm)  [+z into page]")
    ax.set_zlabel("y (mm)  [+y up]")

# One axis thinner than this makes the intersection a sheet, not a volume.
THIN_MM = 1.5
# Two cylinders count as the same fit if their axes lie this close.
MATE_MM = 3.0

ROW = re.compile(r"^\s*(\d+)\s+(\d+)\s+(coplanar|skin|real)\s+([\d.]+)\s+"
                 r"([\d.]+) x\s+([\d.]+) x\s+([\d.]+)\s+(.+?)\s*$")


def read_parts():
    """solid_id -> NX part number (AI-129641-A_1), from stm_placements.csv.

    The overlap report names solids by their Mu2e PV name (LeadFwall9PV), which
    is what the Geant4 side will call them. Going back to the CAD needs the NX
    part number instead, so carry both: the PV name says what it is, the part
    number says where to find it in the STEP.

    Optional: the CSV may not exist yet if the extractor has not run, and a
    missing part number is worth a "?" on a figure, not a dead script.
    """
    path = os.path.join(OUTDIR, "stm_placements.csv")
    out = {}
    if not os.path.exists(path):
        return out
    try:
        import csv
        with open(path) as fh:
            for row in csv.DictReader(fh):
                out[int(row["solid_id"])] = row.get("part") or "?"
    except (ValueError, KeyError, OSError):
        return {}
    return out


PARTS = {}
MODEL_BOXES = []


def part_of(solid_id):
    return PARTS.get(solid_id, "?")


def label(solid_id, pv):
    """"solid 43  AI-129641-A_1  LeadFwall9PV" -- the three names it has."""
    return "solid %d  %s  %s" % (solid_id, part_of(solid_id), pv)


def read_pairs(path):
    """skin/real rows from check_overlaps.py's report."""
    out = []
    with open(path) as fh:
        for line in fh:
            m = ROW.match(line)
            if not m:
                continue
            a, b, kind, vol, ex, ey, ez, parts = m.groups()
            if kind == "coplanar":
                continue
            out.append({"a": int(a), "b": int(b), "kind": kind,
                        "vol": float(vol), "parts": parts.strip()})
    return out


def tris(shape, tol=0.3):
    """Triangles for Poly3DCollection: Mu2e coordinates, display order.

    Each vertex goes through transform() into Mu2e and then through
    mu2e_axes() into matplotlib's (horizontal, depth, vertical) slots, so the
    geometry lands the same way round as every other plot in this project.
    """
    out = []
    for f in shape.Faces:
        try:
            pts, idx = f.tessellate(tol)
        except Exception:
            continue
        for t in idx:
            out.append([S.mu2e_axes(*mu2e_pt(pts[k])) for k in t])
    return out


def cylinders(solid):
    return [f for f in solid.Faces
            if f.Surface.__class__.__name__ == "Cylinder"]


def mating_cylinders(a, b):
    """The one bore/tube pair these two solids actually share, or None.

    A wall can carry several bores -- solid 43 has them 81 mm apart -- so the
    nearest pair is the fit; taking the first would report the distance between
    unrelated holes as an 'offset'.
    """
    best = None
    for fa in cylinders(a):
        for fb in cylinders(b):
            if abs(abs(fa.Surface.Axis.dot(fb.Surface.Axis)) - 1.0) > 1e-6:
                continue
            ca, cb = fa.Surface.Center, fb.Surface.Center
            ax = fa.Surface.Axis
            d = FreeCAD.Vector(cb.x - ca.x, cb.y - ca.y, cb.z - ca.z)
            perp = d.sub(ax * d.dot(ax)).Length
            if perp > MATE_MM:
                continue
            # Prefer the pair of equal radius: that is the interfering fit.
            score = (abs(fa.Surface.Radius - fb.Surface.Radius), perp)
            if best is None or score < best[0]:
                best = (score, fa, fb, perp)
    if best is None:
        return None
    return best[1], best[2], best[3]


def panel_numbers(ax, lines):
    ax.axis("off")
    ax.text(0.0, 1.0, "\n".join(lines), va="top", ha="left", fontsize=8.2,
            family="monospace", transform=ax.transAxes)


def draw_3d(ax, a, b, inter, title, tol=0.35):
    ax.add_collection3d(Poly3DCollection(tris(a, tol), alpha=0.10,
                                         facecolor="#b8860b",
                                         edgecolor="none"))
    ax.add_collection3d(Poly3DCollection(tris(b, tol), alpha=0.18,
                                         facecolor="#2c3e50",
                                         edgecolor="none"))
    ax.add_collection3d(Poly3DCollection(tris(inter, 0.15), alpha=1.0,
                                         facecolor="#c0392b",
                                         edgecolor="#7b241c", linewidth=0.2))
    alo, ahi = mu2e_bbox(a)
    blo, bhi = mu2e_bbox(b)
    lim = {}
    for k, n in enumerate("xyz"):
        lo, hi = min(alo[k], blo[k]), max(ahi[k], bhi[k])
        pad = max((hi - lo) * 0.04, 1e-6)
        lim[n] = (lo - pad, hi + pad)
    mu2e_3d(ax, lim)
    tidy_3d(ax)
    ax.set_title(title, fontsize=9)


def draw_alone(ax, inter, title):
    """The intersection solid on its own, same convention."""
    ax.add_collection3d(Poly3DCollection(tris(inter, 0.15), alpha=0.95,
                                         facecolor="#c0392b",
                                         edgecolor="#7b241c", linewidth=0.25))
    lo, hi = mu2e_bbox(inter)
    lim = {}
    for k, n in enumerate("xyz"):
        pad = max((hi[k] - lo[k]) * 0.06, 1e-6)
        lim[n] = (lo[k] - pad, hi[k] + pad)
    mu2e_3d(ax, lim)
    tidy_3d(ax)
    ax.set_title(title, fontsize=9)


def add_context(fig, pos, solid_ids):
    """Inset: where these two parts sit in the whole model.

    Reuses the same Mu2e-frame boxes plot_all_components draws from, so this
    inset and the one beside a part plot are the same picture of the same
    house. Building it from the CSVs also avoids tessellating 224 solids per
    figure, which would dominate the run.
    """
    boxes = MODEL_BOXES
    ax = fig.add_subplot(*pos, projection="3d")
    ax.set_facecolor("none")
    if not boxes:
        ax.axis("off")
        return ax

    def edges(b):
        hx, hy, hz = b["dx"] / 2.0, b["dy"] / 2.0, b["dz"] / 2.0
        x, y, z = b["x"], b["y"], b["z"]
        c = [(x - hx, y - hy, z - hz), (x + hx, y - hy, z - hz),
             (x + hx, y + hy, z - hz), (x - hx, y + hy, z - hz),
             (x - hx, y - hy, z + hz), (x + hx, y - hy, z + hz),
             (x + hx, y + hy, z + hz), (x - hx, y + hy, z + hz)]
        for u, v in [(0, 1), (1, 2), (2, 3), (3, 0), (4, 5), (5, 6),
                     (6, 7), (7, 4), (0, 4), (1, 5), (2, 6), (3, 7)]:
            pa, pb = S.mu2e_axes(*c[u]), S.mu2e_axes(*c[v])
            yield ([pa[0], pb[0]], [pa[1], pb[1]], [pa[2], pb[2]])

    # model_boxes() carries solid_id as a STRING (it comes straight from the
    # CSV), while the pair ids are ints parsed from the overlap report. Compare
    # as ints on both sides: matching the raw values silently selected nothing,
    # so the highlight loop never ran and the inset drew only the faint house.
    def sid(b):
        try:
            return int(b["solid_id"])
        except (TypeError, ValueError):
            return None

    want = {int(s) for s in solid_ids}
    for b in boxes:
        if sid(b) in want:
            continue
        for xs, ys, zs in edges(b):
            ax.plot(xs, ys, zs, color="#b0b7bd", lw=0.25, alpha=0.45)

    # The pair itself: SOLID faces, not wireframe. A filled box survives being
    # overdrawn by hairlines in a way an outline does not, and at inset scale
    # a 50 mm brick's outline is only a few pixels across either way.
    for b in boxes:
        if sid(b) not in want:
            continue
        col = "#b8860b" if sid(b) == int(solid_ids[0]) else "#2c3e50"
        hx, hy, hz = b["dx"] / 2.0, b["dy"] / 2.0, b["dz"] / 2.0
        x, y, z = b["x"], b["y"], b["z"]
        c = [(x - hx, y - hy, z - hz), (x + hx, y - hy, z - hz),
             (x + hx, y + hy, z - hz), (x - hx, y + hy, z - hz),
             (x - hx, y - hy, z + hz), (x + hx, y - hy, z + hz),
             (x + hx, y + hy, z + hz), (x - hx, y + hy, z + hz)]
        quads = [(0, 1, 2, 3), (4, 5, 6, 7), (0, 1, 5, 4),
                 (3, 2, 6, 7), (0, 3, 7, 4), (1, 2, 6, 5)]
        faces = [[S.mu2e_axes(*c[i]) for i in q] for q in quads]
        ax.add_collection3d(Poly3DCollection(
            faces, facecolor=col, edgecolor=col, linewidth=0.8, alpha=0.95))
        # And a leader from above, so the eye finds the part before it has to
        # pick a small box out of a large house.
        p = S.mu2e_axes(x, y, z)
        top = max(bb["y"] + bb["dy"] / 2.0 for bb in boxes)
        pt = S.mu2e_axes(x, top, z)
        ax.plot([p[0], pt[0]], [p[1], pt[1]], [p[2], pt[2]],
                color=col, lw=1.0, alpha=0.75, ls=":")
        ax.plot([pt[0]], [pt[1]], [pt[2]], marker="v", color=col, ms=6,
                mew=0.6, mec="white")

    lim = {}
    for k in "xyz":
        lo = min(b[k] - b["d" + k] / 2.0 for b in boxes)
        hi = max(b[k] + b["d" + k] / 2.0 for b in boxes)
        lim[k] = (lo, hi)
    mu2e_3d(ax, lim)
    for a in (ax.xaxis, ax.yaxis, ax.zaxis):
        a.set_ticklabels([])
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.set_zlabel("")
    ax.grid(False)
    # The house is 224 faint boxes and the pair is two of them, so at inset
    # size the highlight vanishes unless it is drawn much heavier than the
    # context. Zooming instead would lose the "where in the model" answer the
    # panel exists to give.
    ax.set_title("where they are in the model\n"
                 "gold = solid %d, slate = solid %d" % solid_ids,
                 fontsize=7.5, color="#444444")
    return ax


def figure_tubular(r, a, b, inter, fit):
    """Tube in a bore: the defect is a sub-mm axis offset."""
    fa, fb, perp = fit
    ilo, ihi = mu2e_bbox(inter)
    # Cylinder centres in MU2E, not CAD. The bore axis here is world z in both
    # frames (the transform negates z, which flips the direction but keeps the
    # axis), so the section is still the x-y plane -- but every coordinate on
    # it moves by hundreds of mm, and the offset's x-component changes sign
    # with the flip. Taking these from CAD is what made the earlier figures
    # unreadable next to the model plots.
    ca = mu2e_pt(fa.Surface.Center)
    cb = mu2e_pt(fb.Surface.Center)
    ra, rb = fa.Surface.Radius, fb.Surface.Radius
    dx, dy = cb[0] - ca[0], cb[1] - ca[1]
    off = math.hypot(dx, dy)

    class _P(object):
        def __init__(self, t):
            self.x, self.y, self.z = t
    ca, cb = _P(ca), _P(cb)

    pv_a, pv_b = r["parts"].split(" / ")[0], r["parts"].split(" / ")[-1]
    fig = plt.figure(figsize=(18.5, 8.8))
    fig.suptitle("%s overlap  \u2014  %.1f mm\u00b3\n%s      vs      %s"
                 % (r["kind"].upper(), inter.Volume,
                    label(r["a"], pv_a), label(r["b"], pv_b)),
                 fontsize=12, y=0.985)

    ax1 = fig.add_subplot(2, 4, 1)
    ax1.add_patch(Circle((ca.x, ca.y), ra, fill=False, ec="#b8860b", lw=2.0,
                         label="%s  r=%.3f" % (part_of(r["a"]), ra)))
    ax1.add_patch(Circle((cb.x, cb.y), rb, fill=False, ec="#2c3e50", lw=2.0,
                         label="%s  r=%.3f" % (part_of(r["b"]), rb)))
    ax1.plot([ca.x], [ca.y], "+", color="#b8860b", ms=11, mew=2)
    ax1.plot([cb.x], [cb.y], "+", color="#2c3e50", ms=11, mew=2)
    ax1.set_aspect("equal")
    lim = max(ra, rb) * 1.3
    ax1.set_xlim(ca.x - lim, ca.x + lim)
    ax1.set_ylim(ca.y - lim, ca.y + lim)
    ax1.set_xlabel("Mu2e x (mm)")
    ax1.set_ylabel("Mu2e y (mm)")
    ax1.set_title("section across the axis, TRUE SCALE\n"
                  "at this scale the fit looks concentric", fontsize=9)
    ax1.legend(fontsize=7, loc="lower right")
    ax1.grid(alpha=0.3, lw=0.5)

    ax2 = fig.add_subplot(2, 4, 2)
    ax2.add_patch(Circle((ca.x, ca.y), ra, fill=False, ec="#b8860b", lw=2.0))
    ax2.add_patch(Circle((cb.x, cb.y), rb, fill=False, ec="#2c3e50", lw=2.0))
    ax2.plot([ca.x], [ca.y], "+", color="#b8860b", ms=14, mew=2.5)
    ax2.plot([cb.x], [cb.y], "+", color="#2c3e50", ms=14, mew=2.5)
    if off > 1e-9:
        ax2.annotate("", xy=(cb.x, cb.y), xytext=(ca.x, ca.y),
                     arrowprops=dict(arrowstyle="<->", color="#c0392b",
                                     lw=1.8))
        ax2.text(ca.x + dx / 2.0, ca.y + dy / 2.0,
                 "   %.3f mm\n   (dx %+.3f, dy %+.3f)" % (off, dx, dy),
                 color="#c0392b", fontsize=9, va="center")
    z = max(off * 2.2, 0.6)
    ax2.set_aspect("equal")
    ax2.set_xlim(ca.x - z, ca.x + z)
    ax2.set_ylim(ca.y - z, ca.y + z)
    ax2.set_xlabel("Mu2e x (mm)")
    ax2.set_ylabel("Mu2e y (mm)")
    ax2.set_title("zoom on the two axes (\u00b1%.2f mm)\n"
                  "the defect: a %.3f mm eccentricity" % (z, off), fontsize=9)
    ax2.grid(alpha=0.3, lw=0.5)

    ax3 = fig.add_subplot(2, 4, 3)
    th = [i * math.pi / 400.0 for i in range(801)]
    ax3.fill([ca.x + ra * math.cos(t) for t in th],
             [ca.y + ra * math.sin(t) for t in th],
             color="#b8860b", alpha=0.16)
    ax3.fill([cb.x + rb * math.cos(t) for t in th],
             [cb.y + rb * math.sin(t) for t in th],
             color="#2c3e50", alpha=0.16)
    sx, sy = [], []
    for t in th:
        px, py = cb.x + rb * math.cos(t), cb.y + rb * math.sin(t)
        if math.hypot(px - ca.x, py - ca.y) > ra:
            sx.append(px)
            sy.append(py)
    if sx:
        ax3.plot(sx, sy, color="#c0392b", lw=4, solid_capstyle="round",
                 label="%s outside %s" % (part_of(r["b"]), part_of(r["a"])))
        ax3.legend(fontsize=7, loc="lower right")
    ax3.plot([ca.x + ra * math.cos(t) for t in th],
             [ca.y + ra * math.sin(t) for t in th], color="#b8860b", lw=1.5)
    ax3.plot([cb.x + rb * math.cos(t) for t in th],
             [cb.y + rb * math.sin(t) for t in th], color="#2c3e50", lw=1.5)
    ax3.set_aspect("equal")
    ax3.set_xlim(ca.x - lim, ca.x + lim)
    ax3.set_ylim(ca.y - lim, ca.y + lim)
    ax3.set_xlabel("Mu2e x (mm)")
    ax3.set_ylabel("Mu2e y (mm)")
    ax3.set_title("where one breaches the other\n"
                  "a crescent, deepest along the offset", fontsize=9)
    ax3.grid(alpha=0.3, lw=0.5)

    add_context(fig, (2, 4, 4), (r["a"], r["b"]))

    draw_3d(fig.add_subplot(2, 4, 5, projection="3d"), a, b, inter,
            "solid %d (gold) + solid %d (slate)\n+ overlap (red), in context"
            % (r["a"], r["b"]))

    draw_alone(fig.add_subplot(2, 4, 6, projection="3d"), inter,
               "the overlap solid alone\n%.1f mm\u00b3,  %.2f x %.2f x %.2f"
               % (inter.Volume, ihi[0] - ilo[0], ihi[1] - ilo[1],
                  ihi[2] - ilo[2]))

    panel_numbers(fig.add_subplot(2, 4, (7, 8)), [
        "TUBE / BORE FIT",
        "  solid %-4d %-16s %-22s r %.3f"
        % (r["a"], part_of(r["a"]), pv_a, ra),
        "  solid %-4d %-16s %-22s r %.3f"
        % (r["b"], part_of(r["b"]), pv_b, rb),
        "",
        "  axis A   (%9.3f, %9.3f)  Mu2e" % (ca.x, ca.y),
        "  axis B   (%9.3f, %9.3f)  Mu2e" % (cb.x, cb.y),
        "  offset   dx %+.3f  dy %+.3f  |d| %.3f mm" % (dx, dy, off),
        "  radial clearance  %.3f mm" % (ra - rb),
        "",
        "  intersection %.1f mm\u00b3, %d piece(s)" % (inter.Volume,
                                                       len(inter.Solids)),
        "  extent       %.2f x %.2f x %.2f mm" % (ihi[0] - ilo[0],
                                                  ihi[1] - ilo[1],
                                                  ihi[2] - ilo[2]),
        "  Mu2e x %.1f..%.1f  y %.1f..%.1f  z %.1f..%.1f"
        % (ilo[0], ihi[0], ilo[1], ihi[1], ilo[2], ihi[2]),
        "",
        "READING IT",
        "  Equal radii mean zero clearance, so ANY offset",
        "  interferes. The overlap is a thin crescent running",
        "  the full depth -- large in area, %.3f mm deep." % off,
        "  Invisible in a CAD view at true scale.",
    ])
    fig.tight_layout(rect=(0, 0, 1, 0.925), w_pad=2.6, h_pad=2.2)
    return fig


def figure_planar(r, a, b, inter):
    """Two nominally-coincident faces that actually interpenetrate.

    All coordinates are Mu2e, like every other figure here: the 2D panels name
    the Mu2e axis they section, so a span read off this figure can be looked up
    directly in stm_placements.csv or in the model plots.
    """
    lo, hi = mu2e_bbox(inter)
    ext = [hi[t] - lo[t] for t in range(3)]
    # Which axis is thin is a property of the sheet, not of the frame, but
    # derive it AFTER transforming rather than assuming the CAD answer carries
    # over: the transform flips two axes, and an assumption here would silently
    # section the wrong plane.
    k = min(range(3), key=lambda t: ext[t])
    names = "xyz"
    wide = [t for t in range(3) if t != k]

    pv_a, pv_b = r["parts"].split(" / ")[0], r["parts"].split(" / ")[-1]
    fig = plt.figure(figsize=(18.5, 8.8))
    fig.suptitle("%s overlap  \u2014  %.1f mm\u00b3\n%s      vs      %s"
                 % (r["kind"].upper(), inter.Volume,
                    label(r["a"], pv_a), label(r["b"], pv_b)),
                 fontsize=12, y=0.985)

    # Face-on: the sheet's own plane, showing how much area is shared.
    ax1 = fig.add_subplot(2, 4, 1)
    ax1.add_patch(plt.Rectangle((lo[wide[0]], lo[wide[1]]),
                                ext[wide[0]], ext[wide[1]],
                                facecolor="#c0392b", alpha=0.45,
                                edgecolor="#7b241c", lw=1.5))
    for s, col, tag in ((a, "#b8860b", r["a"]), (b, "#2c3e50", r["b"])):
        slo, shi = mu2e_bbox(s)
        sex = [shi[t] - slo[t] for t in range(3)]
        ax1.add_patch(plt.Rectangle((slo[wide[0]], slo[wide[1]]),
                                    sex[wide[0]], sex[wide[1]], fill=False,
                                    edgecolor=col, lw=1.6,
                                    label=part_of(tag)))
    ax1.set_aspect("equal")
    ax1.relim()
    ax1.autoscale()
    ax1.set_xlabel("Mu2e %s (mm)" % names[wide[0]])
    ax1.set_ylabel("Mu2e %s (mm)" % names[wide[1]])
    ax1.set_title("face-on (%s%s): the shared area\n%.1f x %.1f mm of contact"
                  % (names[wide[0]], names[wide[1]], ext[wide[0]],
                     ext[wide[1]]), fontsize=9)
    ax1.legend(fontsize=7, loc="upper right")
    ax1.grid(alpha=0.3, lw=0.5)

    # Edge-on, zoomed: the whole defect is the depth along the thin axis.
    ax2 = fig.add_subplot(2, 4, 2)
    for s, col, tag in ((a, "#b8860b", r["a"]), (b, "#2c3e50", r["b"])):
        slo, shi = mu2e_bbox(s)
        sex = [shi[t] - slo[t] for t in range(3)]
        ax2.add_patch(plt.Rectangle((slo[wide[0]], slo[k]), sex[wide[0]],
                                    sex[k], facecolor=col, alpha=0.22,
                                    edgecolor=col, lw=1.6,
                                    label=part_of(tag)))
    ax2.add_patch(plt.Rectangle((lo[wide[0]], lo[k]), ext[wide[0]], ext[k],
                                facecolor="#c0392b", alpha=0.85,
                                edgecolor="#7b241c", lw=1.2))
    mid = lo[wide[0]] + ext[wide[0]] / 2.0
    ax2.annotate("", xy=(mid, lo[k]), xytext=(mid, lo[k] + ext[k]),
                 arrowprops=dict(arrowstyle="<->", color="#c0392b", lw=1.8))
    ax2.text(mid, lo[k] + ext[k] / 2.0, "  %.4f mm" % ext[k],
             color="#c0392b", fontsize=10, va="center")
    pad = max(ext[k] * 6.0, 2.0)
    ax2.set_xlim(lo[wide[0]] - ext[wide[0]] * 0.05,
                 lo[wide[0]] + ext[wide[0]] * 1.05)
    ax2.set_ylim(lo[k] - pad, lo[k] + ext[k] + pad)
    ax2.set_xlabel("Mu2e %s (mm)" % names[wide[0]])
    ax2.set_ylabel("Mu2e %s (mm)" % names[k])
    ax2.set_title("edge-on (%s%s), %s stretched\n"
                  "the two solids interpenetrate by %.4f mm"
                  % (names[wide[0]], names[k], names[k], ext[k]), fontsize=9)
    ax2.grid(alpha=0.3, lw=0.5)

    draw_alone(fig.add_subplot(2, 4, 3, projection="3d"), inter,
               "the overlap sheet alone\n(its own aspect: a sheet, not a "
               "block)")

    add_context(fig, (2, 4, 4), (r["a"], r["b"]))

    draw_3d(fig.add_subplot(2, 4, 5, projection="3d"), a, b, inter,
            "solid %d (gold) + solid %d (slate)\n+ overlap (red), in context"
            % (r["a"], r["b"]))

    panel_numbers(fig.add_subplot(2, 4, (6, 8)), [
        "TWO FACES THAT SHOULD MEET, NOT MERGE",
        "  solid %-4d %-16s %s" % (r["a"], part_of(r["a"]), pv_a),
        "  solid %-4d %-16s %s" % (r["b"], part_of(r["b"]), pv_b),
        "",
        "  interpenetration  %.4f mm along Mu2e %s" % (ext[k], names[k]),
        "  shared area       %.1f x %.1f mm" % (ext[wide[0]], ext[wide[1]]),
        "  intersection      %.1f mm\u00b3, %d piece(s)"
        % (inter.Volume, len(inter.Solids)),
        "  Mu2e %s span      %.3f .. %.3f" % (names[k], lo[k], hi[k]),
        "  Mu2e x %.1f..%.1f  y %.1f..%.1f  z %.1f..%.1f"
        % (lo[0], hi[0], lo[1], hi[1], lo[2], hi[2]),
        "",
        "READING IT",
        "  The volume is large only because the area is large.",
        "  The defect is the %.4f mm depth: two faces meant to" % ext[k],
        "  sit against each other are that far inside each other.",
        "  Geant4 will warn; nothing is grossly misplaced.",
    ])
    fig.tight_layout(rect=(0, 0, 1, 0.925), w_pad=2.6, h_pad=2.2)
    return fig


def figure_generic(r, a, b, inter):
    """Neither a sheet nor a tube fit: show it without assuming structure."""
    lo, hi = mu2e_bbox(inter)
    pv_a, pv_b = r["parts"].split(" / ")[0], r["parts"].split(" / ")[-1]
    fig = plt.figure(figsize=(18.5, 5.8))
    fig.suptitle("%s overlap  \u2014  %.1f mm\u00b3\n%s      vs      %s"
                 % (r["kind"].upper(), inter.Volume,
                    label(r["a"], pv_a), label(r["b"], pv_b)),
                 fontsize=12, y=0.975)

    draw_3d(fig.add_subplot(1, 4, 1, projection="3d"), a, b, inter,
            "solid %d (gold) + solid %d (slate)\n+ overlap (red), in context"
            % (r["a"], r["b"]))

    draw_alone(fig.add_subplot(1, 4, 2, projection="3d"), inter,
               "the overlap solid alone\n%.1f mm\u00b3, %.2f x %.2f x %.2f"
               % (inter.Volume, hi[0] - lo[0], hi[1] - lo[1], hi[2] - lo[2]))

    add_context(fig, (1, 4, 3), (r["a"], r["b"]))

    panel_numbers(fig.add_subplot(1, 4, 4), [
        "OVERLAP",
        "  solid %-4d %-16s %s" % (r["a"], part_of(r["a"]), pv_a),
        "  solid %-4d %-16s %s" % (r["b"], part_of(r["b"]), pv_b),
        "",
        "  intersection %.1f mm\u00b3, %d piece(s)" % (inter.Volume,
                                                       len(inter.Solids)),
        "  extent       %.2f x %.2f x %.2f mm" % (hi[0] - lo[0],
                                                  hi[1] - lo[1],
                                                  hi[2] - lo[2]),
        "  Mu2e x span  %.3f .. %.3f" % (lo[0], hi[0]),
        "  Mu2e y span  %.3f .. %.3f" % (lo[1], hi[1]),
        "  Mu2e z span  %.3f .. %.3f" % (lo[2], hi[2]),
        "",
        "  Neither a thin sheet nor a tube-in-bore fit, so",
        "  this is drawn without assuming either shape.",
    ])
    fig.tight_layout(rect=(0, 0, 1, 0.90), w_pad=2.6)
    return fig


def safe(text):
    """Filename-safe, but KEEPING the underscore.

    A part number is written AI-129641-A_1 on the CAD, and the whole point of
    putting it in the filename is that the directory can be searched for that
    exact string. Folding the underscore to a dash gave AI-129641-A-1, which
    looks right and matches nothing.
    """
    return re.sub(r"[^A-Za-z0-9_]+", "-", text).strip("-")


def main():
    report = os.path.join(OUTDIR, "overlaps.txt")
    if not os.path.exists(report):
        sys.exit("missing %s -- run scripts/check_overlaps.py --report first"
                 % report)
    pairs = read_pairs(report)
    if not pairs:
        sys.exit("no skin/real pairs in %s" % report)

    os.makedirs(PLOTDIR, exist_ok=True)
    global PARTS, MODEL_BOXES
    PARTS = read_parts()
    # Loaded once: add_context() draws it per figure, and re-reading the CSVs
    # 37 times would cost more than every boolean in this script.
    try:
        import plot_all_components as PAC
        MODEL_BOXES = PAC.model_boxes()
    except Exception as exc:
        print("  WARNING: no model context (%s)" % exc)
        MODEL_BOXES = []
    print("reading %s" % os.path.basename(STP))
    solids = Part.read(STP).Solids
    print("  %d solids; %d skin/real pair(s) to draw; %d part number(s)"
          % (len(solids), len(pairs), len(PARTS)))

    # The filenames carry the NX part numbers too, so the directory itself can
    # be searched by the number written on the CAD -- grep the listing rather
    # than opening figures to find which one covers AI-129641-A_1.
    stale = [f for f in os.listdir(PLOTDIR) if f.endswith(".png")]
    for f in stale:
        os.remove(os.path.join(PLOTDIR, f))
    if stale:
        print("  removed %d figure(s) from the previous run" % len(stale))

    index = []
    # 'real' first, then by volume: the pair that needs attention leads.
    order = sorted(pairs, key=lambda r: (r["kind"] != "real", -r["vol"]))
    for n, r in enumerate(order, 1):
        a, b = solids[r["a"]], solids[r["b"]]
        inter = a.common(b)
        # Mu2e, like everything else here. The layout choice would survive the
        # CAD box (a rotation and reflection preserve extents up to relabelling
        # the axes), but mixing frames in one script is how the earlier figures
        # came to plot CAD numbers under Mu2e-looking labels.
        _ilo, _ihi = mu2e_bbox(inter)
        ext = [_ihi[t] - _ilo[t] for t in range(3)]

        fit = mating_cylinders(a, b)
        if min(ext) < THIN_MM:
            kind, fig = "planar", figure_planar(r, a, b, inter)
        elif fit is not None:
            kind, fig = "tubular", figure_tubular(r, a, b, inter, fit)
        else:
            kind, fig = "generic", figure_generic(r, a, b, inter)

        name = "%s_%03d_%03d_%s_%s.png" % (r["kind"], r["a"], r["b"],
                                           safe(part_of(r["a"])),
                                           safe(part_of(r["b"])))
        path = os.path.join(PLOTDIR, name)
        fig.savefig(path, dpi=120)
        plt.close(fig)
        index.append({"file": name, "a": r["a"], "b": r["b"],
                      "part_a": part_of(r["a"]), "part_b": part_of(r["b"]),
                      "kind": r["kind"], "layout": kind,
                      "vol": round(inter.Volume, 3), "parts": r["parts"]})
        print("  %2d/%d  %-7s %-7s %s" % (n, len(order), r["kind"], kind,
                                          name))

    with open(os.path.join(PLOTDIR, "index.json"), "w") as fh:
        json.dump(index, fh, indent=1)
    print("\n%d figure(s) in %s" % (len(index), PLOTDIR))
    for k in ("real", "skin"):
        got = [i for i in index if i["kind"] == k]
        if got:
            print("  %-5s %d" % (k, len(got)))


if __name__ == "__main__":
    main()
