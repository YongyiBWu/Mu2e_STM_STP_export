# Mu2e STM shield-house STEP export

Turns the NX STEP export of the STM shield house into a deduplicated, classified
geometry inventory in Mu2e coordinates, ready to be written as Geant4 volumes in
Offline.

## Layout

```
STM_STP_files/   the two NX STEP exports (inputs, read-only)
scripts/         the extractor
scripts/investigation/  the analysis scripts that worked out how (see below)
output/          generated CSVs and JSONs
```

## Inputs

| File | Size | What it is |
|---|---|---|
| `F10269585--_1-G4 Shield House Simplified.stp` | 2.2 MB | 224 solids, fillets removed, **no material data** |
| `F10258491--_1-Shield House Square.stp` | 74 MB | 390 solids, the real assembly, **colour-coded** |

The simplified file is the one to export. The full file exists only because it is
the only place material information survives, and it is colour-coded rather than
carrying real material entities.

## The workflow

Four stages, in order. Each is a separate entry point so you can re-run one
without the others, and `run_investigation.py` chains the first two for you.

```
1. investigate   scripts/run_investigation.py          (8 diagnostic stages)
2. extract       scripts/extract_stm_geometry.py       (writes the 4 CSVs)
3. plot          scripts/plot_all_components.py        (46 PNGs into plots/)
4. examine       scripts/check_overlaps.py             (geometry validation)
```

**The short version** — run everything and be prompted before the extraction:

```
"C:/Program Files/FreeCAD 1.1/bin/python.exe" scripts/run_investigation.py
"C:/Program Files/FreeCAD 1.1/bin/python.exe" scripts/plot_all_components.py
"C:/Program Files/FreeCAD 1.1/bin/python.exe" scripts/check_overlaps.py --report
```

`--report` with no argument writes `output/overlaps.txt`; pass a path to put it
elsewhere. Like `placements_by_shape.csv` it is derived data, so it lives in
`output/` beside the CSVs rather than in `plots/`.

### 1. Investigate

`run_investigation.py` runs the eight stages in `scripts/investigation/`, prints
a timing summary, surfaces the tilt verdict, and offers to run the extractor.
Stage 00 must come first on any new geometry — it decides whether a de-tilt is
even valid. See [the investigation section](#scriptsinvestigation) below.

### 2. Extract

`extract_stm_geometry.py` reads both STEP files and writes the four CSVs. It is
the deliverable; everything else exists to check it or to explain it.

### 3. Plot

`plot_all_components.py` writes one PNG per shape to `plots/`, plus a copy of
the model-wide tilt figure as `_model_tilt_state.png`. Each plot marks the
**Geant4 placement point** with a red crosshair on the part view and prints its
Mu2e coordinates, so a plot and its CSV row can be checked against each other.

What that point is, per kind:

| Kind | Origin marked |
|---|---|
| `BOX`, `BOX_HOLE` | the block centre — a `G4Box` is centred on its origin |
| `TUBE` | the tube's axis midpoint — likewise for `G4Tubs` |
| `PRISM` | the **(0,0) of the cap plane at the sweep half-length** |

The prism case is the one to read carefully. `G4ExtrudedSolid` sweeps a 2D
polygon ±`len/2` and its origin is wherever `(0,0)` falls *in that polygon* —
**not** the polygon's centroid. `stm_prisms.csv` stores the outline relative to
the solid's centre, so `(0,0)` is that centre and the `stm_placements.csv` row
is exactly what to hand `G4PVPlacement`. For shape 13 the cap's area centroid
sits 72 mm away from `(0,0)`; marking the centroid would put the label where
Geant4 places nothing.

Single shapes can still be inspected interactively with
`scripts/show_shape.py --shape N` or `--solid N`.

**A bored shape gets both sections**, because neither answers both questions:

- **Perpendicular to the bore axis** — a hole is a circle, so its radius and
  position read directly. Depth is invisible here, and a counterbore looks
  identical to a through-hole.
- **Parallel to the bore axis** — a hole is a rectangle, showing how deep it
  goes and whether it breaks through the far face. Radius is foreshortened.

The perpendicular view is cut across the **bore axis**, not the block's thinnest
extent. Sectioning the thin axis cuts the plane that *contains* the bore axes,
where drawing a hole as a circle invents geometry and makes coaxial bores at
different depths look like separate holes. Shape 8 exposed it: two z-running
bores, counterbored r=12.7 → r=6.985, appeared as four parallel holes.

#### Display convention

Every 3D view is drawn the way Mu2e is conventionally presented:

```
+y up        -x to the right        +z into the page (at an angle)
```

matplotlib's default 3D triad is +z up with +x toward the viewer, and no
`view_init` changes that — its third argument is *always* the vertical axis. So
the scripts pass Mu2e `(x, y, z)` through a `mu2e_axes()` helper that emits
`(x, z, y)`, then invert the horizontal axis. Axis labels name the Mu2e quantity
on each axis and the direction it runs, because an inverted axis is
indistinguishable from a normal one in a still image.

Sanity check: shape 18 (`SteelBwall18PV`, solid 124) sits at z = −524.169, the
most −z solid in the model, so it must appear at the near end of the depth axis.

One trap worth knowing: `output/tilt_state.json` stores **CAD-frame** centres —
stage 00 describes the input geometry, before any transform — so
`plot_tilt_state.py` applies `extract_stm_geometry.transform()` on read. Drawing
that file's coordinates directly puts the whole model ~350 mm off in x and ~370
mm in z, and mirrored. Everything else plots from `stm_placements.csv`, which is
already Mu2e.

Every plot also carries a **context inset** showing where that shape's copies
sit in the whole model — the house faint, this shape's copies in red — and lists
**every** placement with its NX part id, not just the first. Shapes with many
copies (shape 1 has 146) list the first eight and point at the full CSV.

#### Canonical shape frames

Each shape is stored in a canonical frame so that the Geant4 build can start
from the shape and rotate it into place:

| Kind | Canonical frame |
|---|---|
| `BOX` | the two longest extents on x and y, the shortest as height in z |
| `BOX_HOLE` | longest extent on x — bored blocks are built along their bore axis |

`stm_placements.csv` then carries `r11..r33`, the right-handed rotation taking
that canonical frame to Mu2e axes, plus a readable `orientation` column
(`x->y y->-z z->-x`). Bore offsets and axes in `stm_bores.csv` are expressed in
the **same canonical frame**, so they compose with the box directly.

This matters more than it sounds. Dimensions and rotation have to share one
convention or they cannot be composed: deriving them separately meant the
dimensions were ordered by which world axis each own-axis pointed along, while
the rotation was built from the solid's own axis *order*. Applying one to the
other double-counted the orientation. Defining the canonical frame first, and
the rotation as the map out of it, makes them consistent by construction.

#### Placements and orientations

`list_placements.py` lists every shape with all of its copies — position, the
**NX part id** for cross-referencing against the source assembly, bore count,
and the orientation of each copy:

```
"C:/Program Files/FreeCAD 1.1/bin/python.exe" scripts/list_placements.py --shape 8
"C:/Program Files/FreeCAD 1.1/bin/python.exe" scripts/list_placements.py --csv
```

Orientation is **read from `stm_placements.csv`**, not re-derived — the `r11..r33`
columns and the readable `orientation` string (`x->y y->-z z->-x`), which is
enough to write a `G4RotationMatrix`. Copies with no axis-aligned rotation
(wedges, 45 deg parts) carry an empty orientation rather than a guess.

This script used to recompute orientation from each solid's face normals, and
that disagreed with the CSV: the derivation ordered axes as `orthogonal_frame()`
returned them while the dimensions are ordered canonically, so solid 0 read
`x->x y->-y z->z` here against `x->x y->z z->-y` in the CSV. The CSV value is the
one `verify_roundtrip.py` confirms. Because no STEP read is needed any more, the
script is instant; `--no-summary` just omits the closing tally.

### 4. Examine

`check_overlaps.py` tests whether any two placed solids interpenetrate —
something Geant4 punishes with silently wrong physics rather than an error.

It runs two passes: bounding-box screening (24,976 pairs → a few hundred in
~0.1 s), then a real OCC boolean intersection on each survivor. The second pass
is the one that matters, because a bounding-box test alone would be badly wrong
on this model: the 12 prisms fill only part of their envelope and the 15 bored
blocks have material drilled out, so box-on-box reports overlaps between parts
that merely nest.

A raw intersection count is misleading on its own, because CAD parts that are
*meant* to touch still return a non-zero boolean — OCC gives a coincident
surface a tolerance envelope rather than an exact zero. So each hit is
classified:

| Kind | Meaning |
|---|---|
| `coplanar` | the intersection is a sheet — two faces sharing a plane, e.g. stacked bricks |
| `skin` | a thin shell over a coincident **curved** surface, e.g. a tube in a bore cut to the same radius |
| `real` | one solid genuinely occupying another's space |

The `skin` case needs the volume-fraction test rather than a thin-axis test: a
shared cylinder wraps, so its intersection has no thin bounding-box axis. The
tubes are exactly this — `SteelTube25PV` has outer **r = 25.400** and the lead
block bore it passes through is **r = 25.400**, an exact sliding fit.

**Current result: 173 intersecting pairs, 0 real.** Every one is a shared
surface between parts designed to touch. Geant4 may still warn about coincident
surfaces, but no solid occupies another's space.

`--tol` sets the volume in mm³ below which an intersection is ignored (default
1.0, comfortably above OCC's numerical noise). Exit code is 1 only when a
`real` overlap is found — gating on designed contacts would make the check red
on a correct model — so it works as a build gate.

## The main script

```
"C:/Program Files/FreeCAD 1.1/bin/python.exe" scripts/extract_stm_geometry.py
```

`scripts/extract_stm_geometry.py` is the deliverable. It reads both STEP files and
writes three CSVs:

| Output | Rows | Contents |
|---|---|---|
| `output/stm_shapes.csv` | 45 | one row per distinct shape — the block definitions |
| `output/stm_placements.csv` | 224 | one row per solid — which shape goes where |
| `output/stm_bores.csv` | 32 | cylindrical cuts, for the `G4SubtractionSolid` cases |
| `output/stm_prisms.csv` | 58 | cap outlines for the `G4ExtrudedSolid` cases |

### Bore positions are relative

`stm_bores.csv` gives each bore as `dx,dy,dz` **relative to its block's centre**.
That is what `G4SubtractionSolid` needs, and unlike an absolute position it stays
valid for every placement of a shape — the 146 identical bricks share one
definition, so a world position measured from one of them is meaningless for the
other 145. Absolute `x,y,z` is carried alongside for cross-checking against CAD
only; `dx,dy,dz` is the definition.

### Prisms carry no dx/dy/dz

A prism's `dx,dy,dz` and `in_x,in_y,in_z` in `stm_shapes.csv` are **blank on
purpose**. A prism is a swept outline, not a box, and its bounding envelope is
not a buildable size — shape 19 is a wedge with a 2903 mm² cap inside a
76 × 152 × 76 mm envelope, so building it as a `G4Box` would roughly double the
material. Instead those rows carry `cap_area`, `sweep_len` and `cap_verts`, with
the outline itself in `stm_prisms.csv`.

If a future export contains a prism that is *not* an extrusion, it keeps its
envelope dimensions — that really is all that would be known about it — and gets
no `cap_*` values. So a blank `dx` always means "the outline is in
stm_prisms.csv", never "unknown".

### Prisms are exported as real geometry

All 12 prism solids (11 distinct shapes) are **extrusions**: two congruent
parallel caps joined by walls parallel to the sweep axis, which is exactly
`G4ExtrudedSolid`'s model. `stm_prisms.csv` carries one row per outline vertex,
in order:

| Column | Meaning |
|---|---|
| `shape_id`, `seq` | which shape, and vertex order around the cap |
| `u`, `v` | the vertex in the cap's own 2D plane, mm |
| `len`, `area` | sweep length and cap area |
| `axis_*` | sweep direction, Mu2e coordinates |
| `u_*`, `w_*` | the cap plane's 2D basis in 3D, Mu2e coordinates |

#### Where the vertices are, and how to build them

There is no column of 3D vertices: storing them per placement would repeat the
same polygon for every copy of a shape. The outline is stored **once per shape**,
in 2D, relative to the block centre — and the 3D vertices are recovered by
combining it with a row from `stm_placements.csv`:

```
for each (u, v) in the shape's outline, ordered by seq:
    base    = u * u_vec + v * w_vec                  # u_vec = (u_x, u_y, u_z)
    vertex  = placement + base ± (len/2) * axis_vec  # both caps
```

where `placement` is the `x,y,z` of the solid in `stm_placements.csv`. That
yields all `2 × cap_verts` vertices in Mu2e coordinates, relative to
`_STMShieldingRef`.

The `(u,v)` are **offsets from the block centre**, exactly like the `dx,dy,dz`
in `stm_bores.csv` — so the same composition rule applies to both, and a shape's
outline is valid for every one of its placements.

Verified by round-tripping all 12 prism solids: rebuilding from the CSVs and
comparing against the CAD kernel's own vertices agrees to **0.0007 mm worst
case**. Rebuilt in memory without going through the CSVs the agreement is exact
to floating point, so that residual is formatting alone — now dominated by the
placement `x,y,z` at 3 decimal places. The basis columns `axis_*`/`u_*`/`w_*`
are deliberately written to 10 decimal places: at 4 dp they cost 0.02 mm, since
a 1e-4 rad error swings a vertex ~0.025 mm over a 250 mm arm.

The shape itself is exact — all pairwise vertex distances match the kernel — and
every outline is validated against the kernel's cap area by the shoelace
formula, all 11 agreeing to 1e-5 or better.

It must run under **FreeCAD's bundled interpreter**, not a system Python: it is
FreeCAD that puts the OCC kernel DLLs on the search path, and importing `Part`
first fails with `No module named 'Part'`. Tested against FreeCAD 1.1.

Everything else about the method — the spurious 0.0416° tilt in the export, the
CAD→Mu2e transform, the anchor choice and its independent check against
`stm.STM_SSC.offset_Spot`, why dimensions come from projected vertices rather
than bounding boxes, and what the script deliberately does *not* merge — is
documented in that file's module docstring. Read it before changing anything.

Current result:

```
224 solids, 1402 faces -> 45 distinct shapes
BOX 195 | BOX_HOLE 15 | PRISM 12 | TUBE 2
11 shape definitions cover 190 solids (146 lead bricks in one group)
```

## scripts/investigation/

To run all eight in order and then be asked whether to extract:

```
"C:/Program Files/FreeCAD 1.1/bin/python.exe" scripts/run_investigation.py
```

`scripts/run_investigation.py` is a driver, not a stage. It runs each script as a
subprocess (seven of the eight do their work at import time with no `main()`
guard, so importing them would execute them as a side effect), prints a timing
summary, surfaces the tilt verdict, and then prompts before running
`extract_stm_geometry.py`.

| Flag | Effect |
|---|---|
| `--dry-run` | list the stages and stop |
| `--skip-slow` | skip 02, 03, 05, 06 — the four that read the 74 MB assembly |
| `--only 0,4,7` | run just these stages |
| `--yes` | run the extractor at the end without asking |
| `--no-extract` | stop after the investigation |
| `-q` | suppress each stage's own output |

Stage 02 is skipped with a note if stage 01's JSON is missing, since it reads it;
every other stage works from the STEP files alone, so one failure does not stop
the rest. The driver exits non-zero if any stage failed.

These are the exploratory scripts that worked out *how* to do the extraction.
They are **not** part of the pipeline — `extract_stm_geometry.py` does not import
or invoke any of them, and you never need to run them to produce the CSVs. They
are kept because the reasoning behind the material recovery lives here rather
than in the extractor, and because the numbers they print are the evidence that
the method is sound.

Each one writes a JSON into `output/`. Run them the same way as the main script.

| Script | Writes | Purpose |
|---|---|---|
| `00_detect_global_tilt.py` | `tilt.json` *(only when global)* | **Run this first on any new geometry.** Detects a spurious rotation about z and classifies **each solid** as tilted or square, so a piecewise tilt cannot be mistaken for a global one. Writes the angle for the extractor only when a single rotation is correct for the whole model. See below. |
| `01_parse_msb_colours.py` | `full_msb_colours.json` | Parse the STEP styling chain (`STYLED_ITEM → … → COLOUR_RGB`) to get a colour per solid. Text parsing only — FreeCAD drops the presentation layer, so colour cannot be had from the kernel. |
| `02_bind_colours_by_centroid.py` | `full_solids_colour.json` | First attempt at binding those colours to kernel solids, by centroid proximity. **Superseded by 05** — it is ambiguous where solids share a centroid. Kept because its exact/near/miss counts are what showed the method was insufficient. |
| `03_msb_face_signatures.py` | `full_msb_sigs.json` | Compute a coordinate-free signature (face count + surface types) per solid from the STEP graph, and check the classes agree with the kernel. That agreement is what licenses the positional bind in 05. |
| `04_simplified_inventory.py` | `simplified_solids.json` | Inventory the 224 simplified solids with signature and dimensions. Fast; input to 06. |
| `05_rebuild_colour_lut.py` | `full_solid_colour.json`, `dim_colour_lut.json` | The bind that worked: bucket both sides by face count and zip equal-sized buckets. Resolves all 224. Builds the dimension→colour lookup that became `load_colour_lut()` in the extractor. |
| `06_match_simplified_to_full.py` | `simp_to_full.json` | Match each simplified solid to a full-assembly solid by bounding-box overlap — the hop that carries material between files. Kernel-only, independent of 01's text parsing. Prints an honest match-quality histogram. |
| `07_dedup_groups.py` | `dedup_groups.json` | Group solids into distinct shapes ignoring position and rotation. This is what established that the 146 2×4×8 lead bricks are one definition. Fast. |

Read them in numeric order; it is roughly the order the problem was solved in,
and 02 is left in deliberately as the approach that did not work.

### Checking a new geometry for a spurious rotation

The current export carries a **0.041591 deg rotation about z that is not real** —
an artefact of the NX export. Left in place it corrupts every world-aligned
measurement: a 50.800 mm face reads 51.095, and a bounding box overstates the
45 deg solids threefold. A future export may have no tilt, a different one, or
genuinely rotated parts that must *not* be removed.

`00_detect_global_tilt.py` tells those apart. It takes every planar face normal,
drops the ones pointing along z (a rotation about z leaves them untouched), and
reduces the rest to a deviation from the nearest axis. The verdict turns on
whether the off-axis faces **agree on one angle**:

| Histogram | Meaning | Action |
|---|---|---|
| Off-axis faces agree to ~1e-12 deg, **and every solid carries it** | A global export error — parts placed by hand never agree to that precision | Remove it everywhere |
| Off-axis faces agree, but **some solids are already square** | A *piecewise* tilt | Remove it **per solid**; leave the square ones alone |
| Everything already on axis | Clean export | Remove nothing |
| Off-axis faces scatter over degrees | The parts are genuinely rotated | Remove nothing; measure each solid in its own frame |

**This model is the piecewise case**, which is why the detector classifies each
solid before pronouncing on the model. A pooled histogram cannot separate
"every solid tilted by X" from "most solids tilted by X, some square" — both
produce one tight cluster. The per-solid table is printed above the verdict:

```
tilted    208 of 224
square     12 of 224     15, 18, 23, 29, 117, 118, 139, 163, 164, 168, 170, 173
mixed       2 of 224     135, 220   (genuine 45 deg features)
none        2 of 224     157, 215   (tubes: no in-plane normals to judge by)
```

The full assembly agrees solid-for-solid (302 tilted / 50 square, zero
disagreements across matched pairs), so the split is real CAD structure, not an
artefact of simplification. The square solids cluster in three z-planes and
include the SSC on the beamline.

Because of this, **`tilt.json` is not written for this model.** The detector
emits it only when one global angle is genuinely correct; handing a global angle
to the extractor for a piecewise model is precisely what corrupts the square
solids. The extractor falls back to its built-in `TILT_DEG` until the per-solid
de-tilt is implemented.

It then validates by de-rotating and counting how many faces are *still* off
axis. On the current file that drops 898 → 251, and those 251 are the genuine
45 deg features. **A residue that does not drop means the angle is wrong.**

#### The sign of the de-tilt

The rotation applied is **+`TILT_DEG`**, not −`TILT_DEG`. The transform is
rotate → negate x and z → subtract the anchor, and negating x *mirrors* the
vector, which flips the sign of its angle in the xy plane. Pre-rotating by
−tilt therefore lands at −2×tilt instead of 0.

This is worth knowing because the error is nearly invisible: 2×tilt displaces a
unit vector by only 1.05e-6, well inside the 1e-4 axis-alignment tolerance, and
the box paths use bounding-box dimensions rather than the direction vectors. It
surfaces only where a transformed direction reaches a CSV — `ax,ay,az` in
`stm_bores.csv` and `axis_*`/`u_*`/`w_*` in `stm_prisms.csv`, which read
`0.0015` instead of `0`. The cross-check that settles it: across all 1364
planar normals, +tilt leaves 251 off-axis, −tilt leaves 856.

```
"C:/Program Files/FreeCAD 1.1/bin/python.exe" scripts/investigation/00_detect_global_tilt.py [file.stp]
```

Run against both files today:

| File | Verdict |
|---|---|
| Simplified | Global export tilt, 0.0415914427576354 deg, agreement 1.9e-12 |
| Full assembly | No single global tilt — 1975 off-axis faces spread over 26.5 deg |

The full assembly legitimately contains angled parts, so a global de-tilt would
be wrong there. A model can be *both*, and the script says so: if one deviation
dominates and the rest are few, re-run on a subset you believe is axis-aligned.

### How the tilt reaches the extractor

The detector writes `output/tilt.json` and `extract_stm_geometry.py` reads it,
so the angle is never retyped and keeps full double precision:

```json
{"tilt_deg": 0.04159144275763538,
 "source": "F10269585--_1-G4 Shield House Simplified.stp", ...}
```

The extractor's `TILT_DEG` literal is only a **fallback**, used when the file is
absent, unreadable, or not a float. Two conventions matter if you touch this:

- `tilt_deg` is the **magnitude**; the extractor applies the negation itself
  (`_T = math.radians(-TILT_DEG)`). Handing over a signed value would de-tilt
  the wrong way.
- A `tilt.json` whose `source` names a different file than the one being
  extracted is **ignored with a warning**. A tilt belongs to one export, not to
  the project.

Verified: with `tilt.json` present, absent, malformed, wrong-typed, or naming
the wrong source, all three CSVs come out byte-identical — the measured value
and the literal agree to 3e-14 deg.

So on a new export the sequence is: run `00`, read the verdict, and if it
reports a global tilt the extractor picks it up on its next run.

**Runtime warning.** 02, 03, 05 and 06 read the 74 MB full assembly with the OCC
kernel and take minutes each. 04 and 07 touch only the simplified file and finish
in seconds.

### Note on `output/*.json`

The JSONs in `output/` are the artefacts of the investigation stages above, and
every one of them is reproducible: `scripts/run_investigation.py` rewrites the
lot. The extractor neither reads nor writes them — except `tilt.json`, when
stage 00 writes it — so the JSONs and the CSVs cannot drift into disagreement.

Two are worth calling out:

- `tilt_state.json` — the per-solid tilt classification from stage 00. This is
  the input a per-solid de-tilt needs, and what `scripts/plot_tilt_state.py`
  draws. Centres are in the **CAD frame**, before any transform, because the
  file describes the input geometry rather than the output.
- `tilt.json` — written **only** when a single global angle is correct for the
  whole model. Absent for this geometry; see the tilt section above.

## Status

The material column currently holds **CAD colour names** (`granite gray`,
`medium steel`) as placeholders, with a guessed Geant4 material alongside in
`material_hint` and a `material_match` score per solid so a material inferred
from a loose match is visible rather than implied. Replacing those with real
Geant4 materials, and reviewing the suggested volume names, is the remaining
work. Names are suggestions meant to be edited; the side letter in names like
`CopperLwallPV` comes from a position heuristic and is not authoritative.
