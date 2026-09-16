"""
Emit Geant4 construction code for every shape, with a rotation per copy.

    "C:/Program Files/FreeCAD 1.1/bin/python.exe" scripts/emit_g4.py
    "C:/Program Files/FreeCAD 1.1/bin/python.exe" scripts/emit_g4.py --shape 9
    "C:/Program Files/FreeCAD 1.1/bin/python.exe" scripts/emit_g4.py --out output/constructSTM_generated.cc

Writes output/constructSTM_generated.cc by default: a solid definition per
distinct shape, then one G4PVPlacement per copy carrying that copy's own
G4RotationMatrix.

This is a STARTING POINT, not a drop-in. What it gets right is the geometry --
dimensions, bore subtractions, extruded prisms, positions and orientations, all
from the CSVs. What a human still has to decide:

    materials     the CSVs carry a CAD COLOUR NAME (e.g. "granite gray") and a
                  guess (G4_Pb?). Every findOrBuildMaterial call below is that
                  guess and must be reviewed.
    names         suggested, meant to be edited. The side letter in a name like
                  CopperLwallPV comes from a position heuristic.
    mother volume placements are relative to _STMShieldingRef; the enclosing
                  logical volume is left as a parameter.
    the 16 solids whose orientation could not be reduced to an axis-aligned
                  rotation (wedges and 45 deg parts) are emitted with a comment
                  instead of a rotation, because a wrong matrix is worse than
                  an obvious gap.

Why a rotation per copy, not per shape
    stm_shapes.csv stores dx/dy/dz in each shape's OWN frame -- that is what
    lets 146 lead bricks share one definition. The consequence is that the
    copies differ only in placement and rotation, and the rotation is not
    optional: shape 1 has 13 distinct orientations across its 146 copies, and
    shape 0 has 4 across 20. Emitting one rotation per shape would put most of
    the model at right angles to where it belongs.

G4RotationMatrix convention
    The CSV's r11..r33 map the shape's own frame to Mu2e axes (v_world = R *
    v_own). G4PVPlacement takes the rotation of the MOTHER frame relative to
    the daughter, i.e. the inverse, so the emitted matrix is the transpose.
    That inversion is the single easiest thing to get backwards here, so it is
    applied once, in rotation_code(), and noted at each call site.
"""

import argparse
import collections
import csv
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
OUTDIR = os.path.join(ROOT, "output")


def load(name, required=True):
    path = os.path.join(OUTDIR, name)
    if not os.path.exists(path):
        if not required:
            return []
        sys.exit("missing %s -- run the extractor first" % path)
    with open(path) as fh:
        return list(csv.DictReader(fh))


def num(row, key, default=None):
    v = (row.get(key) or "").strip()
    if not v:
        return default
    try:
        return float(v)
    except ValueError:
        return default


def ident(s):
    """A C++-safe identifier from a suggested volume name."""
    out = "".join(c if (c.isalnum() or c == "_") else "_" for c in (s or ""))
    return out or "unnamed"


def material_call(shape):
    """findOrBuildMaterial for the guessed material, with the colour noted."""
    hint = (shape.get("material_hint") or "").strip().rstrip("?")
    colour = (shape.get("material") or "?").strip()
    if not hint or hint == "?":
        return ('  // MATERIAL UNKNOWN -- CAD colour "%s", no guess available\n'
                '  G4Material* mat = nullptr;  // TODO' % colour)
    return ('  // CAD colour "%s" -> %s (a guess: review before use)\n'
            '  G4Material* mat = findOrBuildMaterial("%s");'
            % (colour, hint, hint))


def rotation_code(place, var):
    """A G4RotationMatrix for one copy, or a comment when none is recoverable.

    The CSV rotation takes the shape's own frame to Mu2e. G4PVPlacement wants
    the inverse, so the rows below are the CSV's COLUMNS -- the transpose.
    """
    if not (place.get("r11") or "").strip():
        return ("    // orientation not axis-aligned (%s); rotation omitted\n"
                "    G4RotationMatrix* %s = nullptr;  // TODO: set by hand"
                % (place.get("orientation") or "oblique/skew", var))
    r = [[num(place, "r%d%d" % (i + 1, j + 1), 0.0) for j in range(3)]
         for i in range(3)]
    if r == [[1, 0, 0], [0, 1, 0], [0, 0, 1]]:
        return "    G4RotationMatrix* %s = nullptr;  // identity" % var
    # Transpose: rowY of the G4 matrix is column j of the CSV rotation.
    rows = []
    for j in range(3):
        rows.append("G4ThreeVector(%s, %s, %s)"
                    % tuple("%g" % r[i][j] for i in range(3)))
    return ("    G4RotationMatrix* %s = new G4RotationMatrix();\n"
            "    %s->rotateAxes(%s,\n"
            "    %s            %s,\n"
            "    %s            %s);  // %s"
            % (var, var, rows[0], " " * len(var), rows[1], " " * len(var),
               rows[2], place.get("orientation") or ""))


def solid_code(shape, bores, prism_rows, var):
    """The G4VSolid for one shape."""
    kind = shape["type"]
    name = ident(shape.get("name"))

    if kind == "TUBE":
        return ('  G4VSolid* %s = new G4Tubs("%s",\n'
                '      %.4f*CLHEP::mm, %.4f*CLHEP::mm, %.4f*CLHEP::mm/2.0,\n'
                '      0.0, CLHEP::twopi);'
                % (var, name, num(shape, "rmin", 0.0), num(shape, "rmax", 0.0),
                   num(shape, "dz", 0.0)))

    if kind in ("PRISM", "OTHER", "PRISM_HOLE") and prism_rows:
        rows = sorted(prism_rows, key=lambda r: int(r["seq"]))
        r0 = rows[0]
        pts = ",\n".join(
            "      G4TwoVector(%.4f*CLHEP::mm, %.4f*CLHEP::mm)"
            % (num(r, "u", 0.0), num(r, "v", 0.0)) for r in rows)
        head = (
            '  // Extruded solid, built the way G4ExtrudedSolid expects: the\n'
            '  // cap polygon lies in x,y and the sweep runs along z, here\n'
            '  // %.4f mm. Nothing about this solid is oblique -- the cap\n'
            '  // basis u=(%.6g,%.6g,%.6g) w=(%.6g,%.6g,%.6g) and the sweep\n'
            '  // axis (%.6g, %.6g, %.6g) form an orthonormal right-handed\n'
            '  // frame, and that frame is the placement rotation below.\n'
            '  std::vector<G4TwoVector> %s_pts = {\n%s };\n'
            '  G4VSolid* %s%s = new G4ExtrudedSolid("%s", %s_pts,\n'
            '      %.4f*CLHEP::mm/2.0,\n'
            '      G4TwoVector(), 1.0, G4TwoVector(), 1.0);'
            % (num(r0, "len", 0.0),
               num(r0, "axis_x", 0.0), num(r0, "axis_y", 0.0),
               num(r0, "axis_z", 0.0),
               num(r0, "u_x", 0.0), num(r0, "u_y", 0.0), num(r0, "u_z", 0.0),
               num(r0, "w_x", 0.0), num(r0, "w_y", 0.0), num(r0, "w_z", 0.0),
               var, pts, var, "_ext" if bores else "", name, var,
               num(r0, "len", 0.0)))
        if not bores:
            return head
        # A PRISM_HOLE is that outline minus its bores. The offsets in
        # stm_bores.csv are already in this same frame -- cap on x,y, sweep on
        # z -- so they need no rotating here, exactly as for BOX_HOLE.
        out = [head]
        prev = "%s_ext" % var
        ln = num(r0, "len", 0.0)
        for n, b in enumerate(bores):
            r = num(b, "r", 0.0)
            depth = num(b, "depth", 0.0) or ln
            ax = (num(b, "ax", 0.0), num(b, "ay", 0.0), num(b, "az", 1.0))
            out.append(
                '  G4VSolid* %s_bore%d = new G4Tubs("%s_bore%d",\n'
                '      0.0, %.4f*CLHEP::mm, %.4f*CLHEP::mm/2.0*1.01,\n'
                '      0.0, CLHEP::twopi);' % (var, n, name, n, r, depth))
            rot = "nullptr"
            if abs(ax[2]) < 0.9:
                rot = "%s_bore%dRot" % (var, n)
                axis = "CLHEP::HepRotationY(CLHEP::halfpi)" \
                    if abs(ax[0]) > 0.9 else \
                    "CLHEP::HepRotationX(CLHEP::halfpi)"
                out.append('  G4RotationMatrix* %s = new G4RotationMatrix(%s);'
                           % (rot, axis))
            out.append(
                '  G4VSolid* %s_cut%d = new G4SubtractionSolid("%s_cut%d",\n'
                '      %s, %s_bore%d, %s,\n'
                '      G4ThreeVector(%.4f*CLHEP::mm, %.4f*CLHEP::mm, '
                '%.4f*CLHEP::mm));'
                % (var, n, name, n, prev, var, n, rot,
                   num(b, "dx", 0.0), num(b, "dy", 0.0), num(b, "dz", 0.0)))
            prev = "%s_cut%d" % (var, n)
        out.append('  G4VSolid* %s = %s;' % (var, prev))
        return "\n".join(out)

    if kind in ("BOX", "BOX_HOLE"):
        dx = num(shape, "dx", 0.0)
        dy = num(shape, "dy", 0.0)
        dz = num(shape, "dz", 0.0)
        out = ['  G4VSolid* %s_box = new G4Box("%s",\n'
               '      %.4f*CLHEP::mm/2.0, %.4f*CLHEP::mm/2.0, '
               '%.4f*CLHEP::mm/2.0);' % (var, name, dx, dy, dz)]
        if not bores:
            out.append('  G4VSolid* %s = %s_box;' % (var, var))
            return "\n".join(out)
        # Each bore is a G4Tubs subtracted at its offset from the block centre.
        prev = "%s_box" % var
        for n, b in enumerate(bores):
            r = num(b, "r", 0.0)
            depth = num(b, "depth", 0.0)
            ax = (num(b, "ax", 0.0), num(b, "ay", 0.0), num(b, "az", 1.0))
            if depth <= 0:
                depth = abs(ax[0]) * dx + abs(ax[1]) * dy + abs(ax[2]) * dz
            # Over-length so the cut face is clean; half-length for G4Tubs.
            out.append(
                '  G4VSolid* %s_bore%d = new G4Tubs("%s_bore%d",\n'
                '      0.0, %.4f*CLHEP::mm, %.4f*CLHEP::mm/2.0*1.01,\n'
                '      0.0, CLHEP::twopi);' % (var, n, name, n, r, depth))
            rot = "nullptr"
            if abs(ax[2]) < 0.9:
                # A bore along x or y needs the tubs turned onto that axis.
                rot = "%s_bore%dRot" % (var, n)
                axis = "CLHEP::HepRotationY(CLHEP::halfpi)" \
                    if abs(ax[0]) > 0.9 else \
                    "CLHEP::HepRotationX(CLHEP::halfpi)"
                out.append('  G4RotationMatrix* %s = new G4RotationMatrix(%s);'
                           % (rot, axis))
            out.append(
                '  G4VSolid* %s_cut%d = new G4SubtractionSolid("%s_cut%d",\n'
                '      %s, %s_bore%d, %s,\n'
                '      G4ThreeVector(%.4f*CLHEP::mm, %.4f*CLHEP::mm, '
                '%.4f*CLHEP::mm));'
                % (var, n, name, n, prev, var, n, rot,
                   num(b, "dx", 0.0), num(b, "dy", 0.0), num(b, "dz", 0.0)))
            prev = "%s_cut%d" % (var, n)
        out.append('  G4VSolid* %s = %s;' % (var, prev))
        return "\n".join(out)

    return '  // shape %s: kind %s not emitted' % (shape["shape_id"], kind)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--shape", type=int, default=None)
    ap.add_argument("--out", default=os.path.join(
        OUTDIR, "constructSTM_generated.cc"))
    args = ap.parse_args()

    shapes = load("stm_shapes.csv")
    places = load("stm_placements.csv")
    bores = collections.defaultdict(list)
    for b in load("stm_bores.csv"):
        bores[b["solid_id"]].append(b)
    prisms = collections.defaultdict(list)
    for r in load("stm_prisms.csv", required=False):
        prisms[r["shape_id"]].append(r)

    by_shape = collections.defaultdict(list)
    for p in places:
        by_shape[p["shape_id"]].append(p)

    if "orientation" not in (places[0] if places else {}):
        sys.exit("stm_placements.csv has no rotation columns -- re-run the "
                 "extractor so copies can be oriented")

    L = []
    L.append("// Generated by scripts/emit_g4.py from the STM STEP export.")
    L.append("// Geometry is derived; MATERIALS AND NAMES ARE GUESSES.")
    L.append("// Positions are relative to _STMShieldingRef, in mm.")
    L.append("//")
    L.append("// Rotations: stm_placements.csv stores R mapping the shape's")
    L.append("// own frame to Mu2e axes. G4PVPlacement wants the inverse, so")
    L.append("// the matrices below are transposed -- see rotation_code().")
    L.append("")
    L.append("#include \"Geant4/G4Box.hh\"")
    L.append("#include \"Geant4/G4Tubs.hh\"")
    L.append("#include \"Geant4/G4ExtrudedSolid.hh\"")
    L.append("#include \"Geant4/G4SubtractionSolid.hh\"")
    L.append("#include \"Geant4/G4PVPlacement.hh\"")
    L.append("#include \"Geant4/G4RotationMatrix.hh\"")
    L.append("")

    n_shapes = 0
    n_copies = 0
    n_norot = 0
    for sh in shapes:
        sid = sh["shape_id"]
        if args.shape is not None and sid != str(args.shape):
            continue
        copies = by_shape.get(sid, [])
        if not copies:
            continue
        n_shapes += 1
        name = ident(sh.get("name"))
        var = "s%s" % sid

        orients = collections.Counter(c.get("orientation") or "(none)"
                                      for c in copies)
        L.append("// " + "-" * 68)
        L.append("// shape %s  %s  %s   x%d"
                 % (sid, sh["type"], name, len(copies)))
        if len(orients) > 1:
            L.append("// %d distinct orientations among its copies:"
                     % len(orients))
            for o, n in orients.most_common():
                L.append("//    %-28s x%d" % (o, n))
        L.append("{")
        L.append(solid_code(sh, bores.get(copies[0]["solid_id"], []),
                            prisms.get(sid, []), var))
        L.append(material_call(sh))
        L.append('  G4LogicalVolume* %sLV = new G4LogicalVolume(%s, mat, "%s");'
                 % (var, var, name))
        L.append("")
        for copy_no, c in enumerate(copies):
            n_copies += 1
            rvar = "%s_rot%s" % (var, c["solid_id"])
            L.append("  {  // solid %s, NX part %s"
                     % (c["solid_id"], c.get("part") or "?"))
            L.append(rotation_code(c, rvar))
            if not (c.get("r11") or "").strip():
                n_norot += 1
            # Name and copy number come from a PER-SHAPE counter, not the
            # solid_id. The solid_id is an index into the STEP file, so using
            # it leaves gaps -- shape 0's twenty copies would be numbered 0, 7,
            # 16, 32, 34, ... which reads as missing volumes and makes the copy
            # number useless for indexing. Geant4 expects 0..n-1 within a
            # logical volume; the NX part id in the comment above is what ties
            # a placement back to the source.
            L.append('    new G4PVPlacement(%s,\n'
                     '        G4ThreeVector(%.3f*CLHEP::mm, %.3f*CLHEP::mm, '
                     '%.3f*CLHEP::mm),\n'
                     '        %sLV, "%s_%d", motherLV, false, %d, doSurfaceCheck);'
                     % (rvar, num(c, "x", 0.0), num(c, "y", 0.0),
                        num(c, "z", 0.0), var, name, copy_no, copy_no))
            L.append("  }")
        L.append("}")
        L.append("")

    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    with open(args.out, "w") as fh:
        fh.write("\n".join(L) + "\n")

    print("wrote %s" % args.out)
    print("  %d shape(s), %d placement(s)" % (n_shapes, n_copies))
    if n_norot:
        print("  %d copy/copies have no axis-aligned rotation and are emitted"
              % n_norot)
        print("  with a TODO instead of a matrix -- search for 'set by hand'")
    return 0


if __name__ == "__main__":
    sys.exit(main())
