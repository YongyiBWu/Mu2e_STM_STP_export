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
| `00_detect_global_tilt.py` | `tilt.json` | **Run this first on any new geometry.** Detects whether the export carries a spurious global rotation about z, and writes the measured angle for the extractor to read. See below. |
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
| Off-axis faces agree to ~1e-12 deg | A global export error — parts placed by hand never agree to that precision | Remove it |
| Everything already on axis | Clean export | Remove nothing |
| Off-axis faces scatter over degrees | The parts are genuinely rotated | Remove nothing; measure each solid in its own frame |

It then validates by de-rotating and counting how many faces are *still* off
axis. On the current file that drops 898 → 251, and those 251 are the genuine
45 deg features. **A residue that does not drop means the angle is wrong.**

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

The JSONs in `output/` are the artefacts of the investigation above, generated
between 16:11 and 16:38 on 2026-09-13. The CSVs are newer, and the extractor
neither reads nor writes the JSONs, so the two cannot drift into disagreement.
`dim_colour_lut_rebuilt.json` is a deliberate keepsake: a copy of
`dim_colour_lut.json` taken before that file was deleted and regenerated, kept so
the rebuilt LUT could be compared against its predecessor.

If you want `output/` to hold only reproducible pipeline results, the JSONs are
safe to delete — the CSVs rebuild from the STEP files alone.

## Status

The material column currently holds **CAD colour names** (`granite gray`,
`medium steel`) as placeholders, with a guessed Geant4 material alongside in
`material_hint` and a `material_match` score per solid so a material inferred
from a loose match is visible rather than implied. Replacing those with real
Geant4 materials, and reviewing the suggested volume names, is the remaining
work. Names are suggestions meant to be edited; the side letter in names like
`CopperLwallPV` comes from a position heuristic and is not authoritative.
