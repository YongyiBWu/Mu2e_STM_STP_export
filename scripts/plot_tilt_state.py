"""
Plot which solids carry the export tilt and which do not.

    "C:/Program Files/FreeCAD 1.1/bin/python.exe" scripts/plot_tilt_state.py [--save FILE]

Reads output/tilt_state.json, written by investigation stage 00: for every
solid, whether its in-plane face normals sit at 0 deg (square) or at
TILT_DEG (tilted). Run stage 00 first -- or the whole chain via
scripts/run_investigation.py -- if the file is missing or stale.

Coordinates in that file are CAD-frame centres, before any transform, because
it describes the input geometry rather than the extractor's output.

Why this plot exists
    The NX export's 0.041591 deg tilt is NOT global -- it is per solid. 208 of
    224 solids carry it, 12 do not, 2 are genuinely oblique (45 deg parts) and
    2 have no in-plane normals to judge by. The extractor currently de-tilts
    everything, which corrects the 208 and corrupts the 12.

    Before changing the transform it is worth seeing whether the split has any
    physical meaning -- a sub-assembly, a region of the house, a size class --
    or whether it is scattered, which would point at the exporter rather than
    the design.

Reading it
    Boxes are each solid's axis-aligned bounding box in Mu2e coordinates,
    drawn as wireframe outlines. The two populations get their own panels with
    SHARED axis limits, so their positions are directly comparable; the 12
    square solids would otherwise be invisible among the 208 tilted ones.
    A third panel overlays both.

    These are bounding boxes, not the solids themselves -- for the oblique
    parts the box overstates the shape. Use show_shape.py for real geometry.
"""

import argparse
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
OUTDIR = os.path.join(ROOT, "output")

# Module level, not inside main(): the per-face plotting below runs at module
# scope and needs transform()/transform_dir() too. main() used to import this
# locally, which was enough when only main() transformed anything.
sys.path.insert(0, HERE)
import extract_stm_geometry as E          # noqa: E402

COLOUR = {
    "tilted": "#c0392b",
    "square": "#1a7f37",
    "oblique": "#b7791f",
    "mixed": "#8e44ad",
    "z-parallel": "#7b8794",
}
LABEL = {
    "tilted": "tilted (carries the 0.0416 deg export tilt)",
    "square": "square (already axis-aligned in CAD)",
    "oblique": "genuinely oblique (45 deg features)",
    "mixed": "mixed (straddles the tilt: some faces square, some tilted)",
    # Kept in the table even though stage 00 now resolves both tubes by their
    # mating bores: a future model may hold a z-parallel solid seated in
    # nothing, which stays unclassifiable and must still draw.
    "z-parallel": "z-parallel (rotation about z unobservable on its own faces)",
}


def mu2e_axes(x, y, z):
    """Map Mu2e (x, y, z) onto matplotlib's (horizontal, depth, vertical).

    matplotlib's 3D THIRD argument is always the vertical axis, and no
    view_init can change that -- rotating the camera tilts the scene but leaves
    z upright. So getting +y up means passing y as the third coordinate and z
    as the second, i.e. swapping them at plot time.

    Every call that draws geometry has to go through this, or the picture
    silently mixes conventions.
    """
    return x, z, y


def mu2e_view(ax):
    """Finish a 3D axes drawn through mu2e_axes(): label it and set the view.

        +y up, -x to the right, +z into the page at an angle.

    The data swap in mu2e_axes() puts y on the vertical axis; this inverts the
    horizontal axis so -x runs to the right, and names each axis by the Mu2e
    quantity actually on it rather than by matplotlib's letter for that slot.

    Sanity check: shape 18 (SteelBwall18PV, solid 124) sits at z = -524.169,
    the most -z solid in the model, so it must appear at the near end of the
    depth axis.
    """
    # Raised viewpoint, matching the per-shape plots: looking down from well
    # above +y reads the model far better than a near-horizontal view, where
    # the shield house's own walls hide everything behind them.
    ax.view_init(elev=32, azim=-68)
    ax.set_xlim(ax.get_xlim()[::-1])          # -x to the right
    ax.set_xlabel("x (mm)  [-x right]", fontsize=8)
    ax.set_ylabel("z (mm)  [+z into page]", fontsize=8)
    ax.set_zlabel("y (mm)  [+y up]", fontsize=8)


def box_edges(x, y, z, dx, dy, dz):
    """The 12 edges of an axis-aligned box, as (xs, ys, zs) triples."""
    hx, hy, hz = dx / 2.0, dy / 2.0, dz / 2.0
    c = [(x - hx, y - hy, z - hz), (x + hx, y - hy, z - hz),
         (x + hx, y + hy, z - hz), (x - hx, y + hy, z - hz),
         (x - hx, y - hy, z + hz), (x + hx, y - hy, z + hz),
         (x + hx, y + hy, z + hz), (x - hx, y + hy, z + hz)]
    idx = [(0, 1), (1, 2), (2, 3), (3, 0),
           (4, 5), (5, 6), (6, 7), (7, 4),
           (0, 4), (1, 5), (2, 6), (3, 7)]
    # Emit in matplotlib's axis order via mu2e_axes, so the vertical axis
    # carries Mu2e y.
    out = []
    for a, b in idx:
        pa = mu2e_axes(*c[a])
        pb = mu2e_axes(*c[b])
        out.append(([pa[0], pb[0]], [pa[1], pb[1]], [pa[2], pb[2]]))
    return out


def draw_set(ax, rows, states, limits, title, lw=0.6, alpha=0.55):
    for r in rows:
        if r["state"] not in states:
            continue
        col = COLOUR.get(r["state"], "#000000")
        for xs, ys, zs in box_edges(r["x"], r["y"], r["z"],
                                    r["dx"], r["dy"], r["dz"]):
            ax.plot(xs, ys, zs, color=col, linewidth=lw, alpha=alpha)
    # limits arrive as Mu2e (x, y, z); reorder to match the plotted axes.
    lx, ly, lz = limits
    ax.set_xlim(*lx)
    ax.set_ylim(*lz)
    ax.set_zlim(*ly)
    try:
        ax.set_box_aspect([lx[1] - lx[0], lz[1] - lz[0], ly[1] - ly[0]])
    except (AttributeError, ValueError):
        pass
    mu2e_view(ax)
    ax.tick_params(labelsize=7)
    ax.set_title(title, fontsize=10)


def face_states(solid, tilt_deg, z_cut=0.99, onaxis=1e-6, agree=1e-3):
    """Classify every face of one solid, the way stage 00 classifies solids.

    Returns a list of dicts: index, kind, area, Mu2e normal, deviation from
    the nearest axis in degrees, and the state that deviation implies.

    A face whose normal is within z_cut of +/-z is 'z-parallel': a rotation
    ABOUT z maps it to itself, so it is evidence of nothing either way. Saying
    that outright matters on this plot -- otherwise the two big z-faces look
    like two more 'square' votes.
    """
    import math
    out = []
    for j, f in enumerate(solid.Faces):
        cn = f.Surface.__class__.__name__
        if cn != "Plane":
            ax = f.Surface.Axis
            out.append({"i": j, "kind": cn, "area": f.Area,
                        "n": E.transform_dir(ax.x, ax.y, ax.z),
                        "dev": None, "state": "curved",
                        "centre": None})
            continue
        n = f.normalAt(0, 0)
        n.normalize()
        # Deviation is measured on the RAW CAD normal, like stage 00: the
        # transform removes the very tilt being measured.
        if abs(n.z) > z_cut:
            dev, state = None, "z-parallel"
        else:
            a = math.degrees(math.atan2(n.y, n.x)) % 90.0
            dev = abs(a if a < 45.0 else a - 90.0)
            if dev <= onaxis:
                state = "square"
            elif abs(dev - tilt_deg) <= agree:
                state = "tilted"
            else:
                state = "oblique"
        c = f.CenterOfMass
        out.append({"i": j, "kind": "Plane", "area": f.Area,
                    "n": E.transform_dir(n.x, n.y, n.z),
                    "dev": dev, "state": state,
                    "centre": E.transform(c.x, c.y, c.z)})
    return out


def figure_faces(solid, info, tilt_deg, part, plt):
    """One mixed solid, face by face, with every normal drawn and named."""
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection
    from matplotlib.lines import Line2D

    fig = plt.figure(figsize=(17.5, 7.4))
    fig.suptitle("solid %d  %s  —  MIXED: some faces carry the "
                 "%.6f deg tilt, some are already square"
                 % (info["i"], part, tilt_deg), fontsize=12, y=0.98)

    faces = face_states(solid, tilt_deg)
    planar = [f for f in faces if f["kind"] == "Plane"]
    span = max(solid.BoundBox.XLength, solid.BoundBox.YLength,
               solid.BoundBox.ZLength)
    arrow = span * 0.45

    # --- panel 1: the solid, with a normal drawn out of every planar face ---
    ax = fig.add_subplot(1, 3, 1, projection="3d")
    tri = []
    for f in solid.Faces:
        try:
            pts, idx = f.tessellate(0.3)
        except Exception:
            continue
        for t in idx:
            tri.append([mu2e_axes(*E.transform(pts[k].x, pts[k].y, pts[k].z))
                        for k in t])
    if tri:
        ax.add_collection3d(Poly3DCollection(tri, alpha=0.16,
                                             facecolor="#8e44ad",
                                             edgecolor="none"))
    for f in planar:
        col = COLOUR.get(f["state"], "#7b8794")
        c, n = f["centre"], f["n"]
        tip = (c[0] + n[0] * arrow, c[1] + n[1] * arrow, c[2] + n[2] * arrow)
        pa, pb = mu2e_axes(*c), mu2e_axes(*tip)
        ax.plot([pa[0], pb[0]], [pa[1], pb[1]], [pa[2], pb[2]],
                color=col, lw=2.0, zorder=5)
        ax.plot([pb[0]], [pb[1]], [pb[2]], marker="o", color=col, ms=4,
                mec="white", mew=0.6, zorder=6)
        ax.text(pb[0], pb[1], pb[2], " f%d" % f["i"], color=col, fontsize=7.5)

    pts_all = [mu2e_axes(*E.transform(v.Point.x, v.Point.y, v.Point.z))
               for v in solid.Vertexes]
    for setter, k in ((ax.set_xlim, 0), (ax.set_ylim, 1), (ax.set_zlim, 2)):
        lo = min(p[k] for p in pts_all) - arrow * 1.15
        hi = max(p[k] for p in pts_all) + arrow * 1.15
        setter(lo, hi)
    try:
        ax.set_box_aspect([1, 1, 1])
    except (AttributeError, ValueError):
        pass
    mu2e_view(ax)
    ax.tick_params(labelsize=6.5)
    ax.set_title("every planar face's outward normal\n"
                 "colour = what that face says about the tilt", fontsize=9)

    # --- panel 2: the deviations, where the angle is actually legible -------
    #
    # The arrows in panel 1 CANNOT show the tilt: 0.0416 deg is 0.7 mm of
    # deflection at the end of a 1 m arrow. Panel 1 shows which face is which;
    # this panel is where the angle itself is readable.
    ax2 = fig.add_subplot(1, 3, 2)
    judged = [f for f in planar if f["dev"] is not None]
    ys = list(range(len(judged)))
    for y, f in zip(ys, judged):
        col = COLOUR.get(f["state"], "#7b8794")
        ax2.barh(y, f["dev"], color=col, alpha=0.85, height=0.62)
        ax2.text(f["dev"] + tilt_deg * 0.04, y,
                 "  %.9f deg" % f["dev"], va="center", fontsize=7.5,
                 color=col)
    ax2.axvline(0.0, color="#1a7f37", lw=1.2, ls="--")
    ax2.axvline(tilt_deg, color="#c0392b", lw=1.2, ls="--")
    ax2.text(tilt_deg, len(judged) - 0.35, " export tilt", fontsize=7.5,
             color="#c0392b", va="top")
    ax2.set_yticks(ys)
    ax2.set_yticklabels(["f%d  %.0f mm²" % (f["i"], f["area"])
                         for f in judged], fontsize=7.5)
    ax2.set_xlim(-tilt_deg * 0.12, tilt_deg * 1.75)
    ax2.set_xlabel("deviation of the face normal from the nearest axis (deg)")
    ax2.set_title("the angle itself\n"
                  "0 = square, %.4f = carries the export tilt" % tilt_deg,
                  fontsize=9)
    ax2.grid(axis="x", alpha=0.3, lw=0.5)

    # --- panel 3: the table, and what it means ------------------------------
    ax3 = fig.add_subplot(1, 3, 3)
    ax3.axis("off")
    lines = ["FACE BY FACE   (normals in Mu2e coordinates)", ""]
    lines.append("  face  area mm²   normal (x, y, z)           state")
    for f in faces:
        n = f["n"]
        if f["kind"] != "Plane":
            lines.append("  f%-4d %9.1f  axis (%6.3f %6.3f %6.3f)  %s"
                         % (f["i"], f["area"], n[0], n[1], n[2], f["kind"]))
            continue
        lines.append("  f%-4d %9.1f  (%7.4f %7.4f %7.4f)  %s"
                     % (f["i"], f["area"], n[0], n[1], n[2], f["state"]))
    counts = {}
    for f in planar:
        counts[f["state"]] = counts.get(f["state"], 0) + 1
    lines += ["", "WHY THIS SOLID IS 'MIXED'",
              "  %d face(s) tilted, %d square, %d z-parallel."
              % (counts.get("tilted", 0), counts.get("square", 0),
                 counts.get("z-parallel", 0))]

    # The clearest statement of the defect: two faces on the SAME side of the
    # block disagreeing. Found by grouping planar faces by which world axis
    # they face and looking for a group holding both states.
    import collections
    by_axis = collections.defaultdict(list)
    for f in planar:
        if f["dev"] is None:
            continue
        n = f["n"]
        k = max(range(3), key=lambda t: abs(n[t]))
        by_axis[("xyz"[k], n[k] > 0)].append(f)
    for (axis, pos), group in sorted(by_axis.items()):
        states = {f["state"] for f in group}
        if len(states) > 1:
            lines.append("")
            lines.append("  The %s%s side disagrees with ITSELF:"
                         % ("+" if pos else "-", axis))
            for f in sorted(group, key=lambda f: -f["area"]):
                lines.append("     f%-3d %9.1f mm²  %s"
                             % (f["i"], f["area"], f["state"]))
    lines += ["", "CONSEQUENCE",
              "  A blanket de-tilt rotates every face by the same",
              "  angle, so it squares the tilted faces and knocks",
              "  the square ones off axis by %.6f deg." % tilt_deg,
              "  Whichever way it is applied, part of this solid",
              "  ends up wrong.",
              "",
              "  z-parallel faces are listed for completeness: a",
              "  rotation about z maps them to themselves, so they",
              "  are evidence of nothing either way."]
    ax3.text(0.0, 1.0, "\n".join(lines), va="top", ha="left", fontsize=8.0,
             family="monospace", transform=ax3.transAxes)

    seen = [s for s in ("tilted", "square", "z-parallel", "oblique")
            if any(f["state"] == s for f in planar)]
    fig.legend(handles=[Line2D([0], [0], color=COLOUR.get(s, "#7b8794"), lw=2.5,
                               label=LABEL.get(s, s)) for s in seen],
               loc="lower center", ncol=4, fontsize=8, frameon=False,
               bbox_to_anchor=(0.5, 0.005))
    fig.tight_layout(rect=(0, 0.05, 1, 0.93))
    return fig


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--save", metavar="FILE")
    ap.add_argument("--faces-dir", metavar="DIR",
                    help="where the per-mixed-solid face plots go "
                         "(default: beside --save)")
    ap.add_argument("--no-faces", action="store_true",
                    help="skip the per-mixed-solid face plots")
    args = ap.parse_args()

    path = os.path.join(OUTDIR, "tilt_state.json")
    if not os.path.exists(path):
        sys.exit("missing %s" % path)
    rows = json.load(open(path))

    # tilt_state.json holds CAD-FRAME centres: stage 00 describes the input
    # geometry, before any transform, which is the right thing for a file about
    # what the exporter produced. But this plot is labelled Mu2e coordinates and
    # sits beside plots built from stm_placements.csv, so the centres have to be
    # transformed on the way in. Without this the whole model is drawn shifted
    # by the anchor (~350mm in x, ~370mm in z) and mirrored in x and z.
    #
    # Extents are unsigned lengths; the transform is a rotation by a few
    # thousandths of a degree plus a reflection, so dx/dy/dz carry over as-is.
    # Verified against stm_placements.csv: all 224 solids agree to 0.001 mm.
    for r in rows:
        r["x"], r["y"], r["z"] = E.transform(r["x"], r["y"], r["z"])

    import matplotlib
    if args.save:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    # Shared limits across every panel, so positions compare directly.
    lim = []
    for k in "xyz":
        lo = min(r[k] - r["d" + k] / 2.0 for r in rows)
        hi = max(r[k] + r["d" + k] / 2.0 for r in rows)
        pad = (hi - lo) * 0.04
        lim.append((lo - pad, hi + pad))

    counts = {}
    for r in rows:
        counts[r["state"]] = counts.get(r["state"], 0) + 1

    # NX part numbers, for naming the odd solids below. Optional on purpose:
    # stage 00 can run before the extractor has ever written a CSV, and a
    # missing part number is worth a "?" in a caption, not a dead plot.
    parts = {}
    csv_path = os.path.join(OUTDIR, "stm_placements.csv")
    if os.path.exists(csv_path):
        import csv as _csv
        try:
            with open(csv_path) as fh:
                for row in _csv.DictReader(fh):
                    parts[int(row["solid_id"])] = row.get("part") or "?"
        except (ValueError, KeyError, OSError):
            parts = {}

    fig = plt.figure(figsize=(16.5, 6.2))

    ax1 = fig.add_subplot(1, 3, 1, projection="3d")
    draw_set(ax1, rows, {"tilted"}, lim,
             "TILTED -- %d solids\nde-tilt is correct for these"
             % counts.get("tilted", 0))

    # 'mixed' belongs here, not only in the overlay. A mixed solid has some
    # faces square and some tilted -- solid 220 straddles the split -- so a
    # blanket de-tilt is wrong for part of it whichever way it is applied.
    # Leaving it out of this panel implied the only non-tilted parts were the
    # three square ones.
    ax2 = fig.add_subplot(1, 3, 2, projection="3d")
    draw_set(ax2, rows, {"square", "oblique", "mixed", "z-parallel"}, lim,
             "NOT FULLY TILTED -- %d square, %d mixed, %d oblique, "
             "%d z-parallel\nde-tilt CORRUPTS the square ones and half of "
             "each mixed one"
             % (counts.get("square", 0), counts.get("mixed", 0),
                counts.get("oblique", 0), counts.get("z-parallel", 0)),
             lw=1.3, alpha=0.95)

    ax3 = fig.add_subplot(1, 3, 3, projection="3d")
    draw_set(ax3, rows, set(COLOUR), lim, "both, overlaid", lw=0.5, alpha=0.5)

    handles = [Line2D([0], [0], color=COLOUR[s], lw=2,
                      label="%s  (%d)" % (LABEL[s], counts.get(s, 0)))
               for s in ("tilted", "square", "oblique", "mixed", "z-parallel")
               if counts.get(s)]
    fig.legend(handles=handles, loc="lower center", ncol=4, fontsize=8,
               frameon=False, bbox_to_anchor=(0.5, 0.005))

    # Name the exceptions outright. The tilted set is 219 solids and nobody
    # will read it off a plot, but the handful that are NOT simply tilted are
    # exactly what someone porting this geometry has to look up by hand, so
    # print them -- solid index and NX part number both, since the CAD is
    # navigated by part number and the CSVs by solid id.
    odd = [r for r in rows if r["state"] in ("square", "mixed", "oblique")]
    if odd:
        print("square / mixed / oblique solids (the ones a blanket de-tilt "
              "would damage):")
        for r in sorted(odd, key=lambda r: (r["state"], r["i"])):
            print("  solid %-4d %-8s part %-18s %s"
                  % (r["i"], r["state"], parts.get(r["i"], "?"),
                     r.get("note") or ""))
        caption = "   ".join(
            "%s: %d (%s)" % (r["state"], r["i"], parts.get(r["i"], "?"))
            for r in sorted(odd, key=lambda r: (r["state"], r["i"])))
        fig.text(0.5, 0.055, caption, ha="center", fontsize=7.5,
                 color="#444444")

    fig.suptitle("STM shield house: which solids carry the NX export tilt "
                 "(bounding boxes, Mu2e coordinates)", fontsize=11, y=0.995)
    # Leave room at the top for the suptitle and at the bottom for the legend;
    # without this the suptitle prints straight over the middle panel's title.
    fig.tight_layout(rect=(0, 0.06, 1, 0.93))

    if args.save:
        fig.savefig(args.save, dpi=130)
        print("wrote %s" % args.save)

    # One figure per MIXED solid, face by face.
    #
    # A mixed solid is the only case the summary panels cannot explain: a
    # tilted solid is tilted and a square one is square, but "mixed" is a claim
    # about disagreement WITHIN one part, and the only way to show that is to
    # name the faces. Reading the STEP is a new dependency for this script, so
    # a failure here warns and leaves the main figure standing.
    mixed = [r for r in rows if r["state"] == "mixed"]
    if mixed and not args.no_faces:
        outdir = args.faces_dir or (os.path.dirname(os.path.abspath(args.save))
                                    if args.save else OUTDIR)
        try:
            import Part
            solids = Part.read(E.SIMPLE).Solids
        except Exception as exc:
            print("  WARNING: no per-face plots (%s)" % exc)
            solids = None
        if solids is not None:
            os.makedirs(outdir, exist_ok=True)
            for r in mixed:
                part = parts.get(r["i"], "?")
                fig2 = figure_faces(solids[r["i"]], r, E.TILT_DEG, part, plt)
                name = "tilt_faces_%03d_%s.png" % (
                    r["i"], part.replace("/", "-"))
                path2 = os.path.join(outdir, name)
                fig2.savefig(path2, dpi=130)
                plt.close(fig2)
                print("wrote %s" % path2)

    if not args.save:
        plt.show()


if __name__ == "__main__":
    main()
