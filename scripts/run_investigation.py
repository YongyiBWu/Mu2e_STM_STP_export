"""
Run the eight investigation stages in order, then offer to run the extractor.

    "C:/Program Files/FreeCAD 1.1/bin/python.exe" scripts/run_investigation.py

Options
    --yes           run the extractor at the end without asking
    --no-extract    stop after the investigation, do not ask
    --only 0,4,7    run just these stages (comma-separated, in the given order)
    --skip-slow     skip the stages that read the 74MB full assembly (02,03,05,06)
    --dry-run       list what would run and stop

Why a driver
    The stages are separate scripts on purpose -- each answers one question and
    can be re-run alone while working on it. But run from cold they have an
    order, and one real dependency: 02 reads the JSON that 01 writes. Running
    them by hand means remembering that, and remembering which four take minutes
    because they read the full assembly with the OCC kernel.

    This runs them as SUBPROCESSES rather than importing them. Seven of the
    eight do their work at module level with no main() guard, so importing would
    execute them as a side effect and any failure would take the driver down
    with it. A subprocess also isolates the OCC kernel, which is the part most
    likely to die badly on a malformed STEP file.

What it does NOT do
    It does not re-run a stage whose inputs are missing because an earlier stage
    failed: 02 is skipped if 01 did not produce its JSON. Every other stage
    reads only the STEP files, so a failure there is reported and the run
    continues -- one broken stage should not hide the results of the other six.

The extractor prompt
    The investigation is diagnostic; extract_stm_geometry.py is what actually
    writes the CSVs. They are separated because on a new geometry you want to
    READ the tilt verdict before trusting any extraction. So the default is to
    ask, and to show the tilt verdict in the summary first.
"""

import argparse
import collections
import json
import os
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
INVDIR = os.path.join(HERE, "investigation")
OUTDIR = os.path.join(ROOT, "output")
EXTRACTOR = os.path.join(HERE, "extract_stm_geometry.py")

# (stage number, filename, what it writes, reads the 74MB full assembly?)
STAGES = [
    ("00", "00_detect_global_tilt.py", ["tilt.json"], False),
    ("01", "01_parse_msb_colours.py", ["full_msb_colours.json"], False),
    ("02", "02_bind_colours_by_centroid.py", ["full_solids_colour.json"], True),
    ("03", "03_msb_face_signatures.py", ["full_msb_sigs.json"], True),
    ("04", "04_simplified_inventory.py", ["simplified_solids.json"], False),
    ("05", "05_rebuild_colour_lut.py",
     ["dim_colour_lut.json", "full_solid_colour.json"], True),
    ("06", "06_match_simplified_to_full.py", ["simp_to_full.json"], True),
    ("07", "07_dedup_groups.py", ["dedup_groups.json"], False),
]

# Stage -> the outputs it needs from an earlier stage. Only 02 has one: it reads
# what 01 writes. Everything else works from the STEP files alone.
NEEDS = {"02": ["full_msb_colours.json"]}

Result = collections.namedtuple("Result", "stage name status secs note")


def run_stage(stage, name, verbose):
    """Run one stage as a subprocess. Returns (status, seconds, note)."""
    path = os.path.join(INVDIR, name)
    missing = [f for f in NEEDS.get(stage, [])
               if not os.path.exists(os.path.join(OUTDIR, f))]
    if missing:
        return ("SKIP", 0.0, "needs %s from an earlier stage" % ", ".join(missing))

    t0 = time.time()
    try:
        p = subprocess.Popen([sys.executable, path],
                             stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                             cwd=ROOT, universal_newlines=True)
        out, _ = p.communicate()
        code = p.returncode
    except Exception as exc:                      # pragma: no cover
        return ("FAIL", time.time() - t0, "could not start: %s" % exc)
    secs = time.time() - t0

    if verbose:
        for line in (out or "").rstrip().splitlines():
            print("      | " + line)
    if code != 0:
        tail = [l for l in (out or "").rstrip().splitlines() if l.strip()]
        note = tail[-1][:70] if tail else "exit %d, no output" % code
        return ("FAIL", secs, note)
    return ("OK", secs, "")


def tilt_verdict():
    """Read back what stage 00 concluded, for the summary."""
    try:
        with open(os.path.join(OUTDIR, "tilt.json")) as fh:
            rec = json.load(fh)
    except (IOError, OSError, ValueError):
        return None
    return rec


def main():
    ap = argparse.ArgumentParser(add_help=True, description=__doc__.strip(),
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--yes", action="store_true",
                    help="run the extractor at the end without asking")
    ap.add_argument("--no-extract", action="store_true",
                    help="stop after the investigation")
    ap.add_argument("--only", default=None,
                    help="comma-separated stage numbers, e.g. 0,4,7")
    ap.add_argument("--skip-slow", action="store_true",
                    help="skip stages that read the full assembly")
    ap.add_argument("--dry-run", action="store_true",
                    help="list what would run and stop")
    ap.add_argument("-q", "--quiet", action="store_true",
                    help="suppress each stage's own output")
    args = ap.parse_args()

    stages = STAGES
    if args.only:
        want = {s.strip().zfill(2) for s in args.only.split(",") if s.strip()}
        unknown = want - {s[0] for s in STAGES}
        if unknown:
            ap.error("unknown stage(s): %s" % ", ".join(sorted(unknown)))
        stages = [s for s in STAGES if s[0] in want]
    if args.skip_slow:
        stages = [s for s in stages if not s[3]]
    if not stages:
        print("nothing to run")
        return 0

    print("=" * 62)
    print("STM investigation: %d stage(s)" % len(stages))
    print("  interpreter: %s" % sys.executable)
    print("  output:      %s" % OUTDIR)
    slow = sum(1 for s in stages if s[3])
    if slow:
        print("  %d stage(s) read the 74MB assembly and take minutes each" % slow)
    print("=" * 62)

    if args.dry_run:
        for stage, name, writes, is_slow in stages:
            print("  %s  %-34s -> %-42s%s"
                  % (stage, name, ", ".join(writes), "  [slow]" if is_slow else ""))
        return 0

    results = []
    t_all = time.time()
    for i, (stage, name, writes, is_slow) in enumerate(stages, 1):
        print("\n[%d/%d] %s%s" % (i, len(stages), name, "  (slow)" if is_slow else ""))
        status, secs, note = run_stage(stage, name, not args.quiet)
        results.append(Result(stage, name, status, secs, note))
        mark = {"OK": "ok", "FAIL": "FAILED", "SKIP": "skipped"}[status]
        print("      %s  %.1fs%s" % (mark, secs, "  -- " + note if note else ""))

    print("\n" + "=" * 62)
    print("summary  (%.1fs total)" % (time.time() - t_all))
    print("=" * 62)
    for r in results:
        print("  %s  %-34s %-8s %6.1fs %s"
              % (r.stage, r.name, r.status, r.secs, r.note))

    ok = sum(1 for r in results if r.status == "OK")
    bad = [r for r in results if r.status == "FAIL"]
    print("\n  %d ok, %d failed, %d skipped"
          % (ok, len(bad), sum(1 for r in results if r.status == "SKIP")))

    # The tilt verdict is the one result worth reading before extracting, so
    # surface it rather than leaving it in stage 00's scrollback.
    rec = tilt_verdict()
    print()
    if rec is None:
        print("  tilt: no output/tilt.json -- stage 00 found no global tilt, or")
        print("        did not run. The extractor will use its built-in value.")
    else:
        print("  tilt: %.15g deg measured from %s"
              % (rec.get("tilt_deg", float("nan")), rec.get("source", "?")))
        print("        off-axis faces agree to %.2e deg; de-tilt takes the"
              % rec.get("agreement_deg", float("nan")))
        print("        off-axis count %s -> %s of %s"
              % (rec.get("residue_before", "?"), rec.get("residue_after", "?"),
                 rec.get("planar_faces", "?")))
        print("        the extractor will pick this up automatically.")

    if args.no_extract:
        print("\n--no-extract given; stopping before extract_stm_geometry.py")
        return 1 if bad else 0

    print()
    if bad:
        print("  NOTE: %d stage(s) failed. The extractor reads only the STEP"
              % len(bad))
        print("        files, so it can still run -- but read the failures")
        print("        above first in case they say something about the input.")

    if args.yes:
        go = True
    elif not sys.stdin.isatty():
        print("  not a terminal and --yes not given; stopping. Run the")
        print("  extractor yourself, or pass --yes.")
        return 1 if bad else 0
    else:
        # isatty() is not reliable here: under Git Bash on Windows a redirected
        # stdin can still report as a tty, so the check above lets us through
        # and the read hits EOF. Catching it is what actually makes the
        # non-interactive case behave, so both reads are wrapped.
        prompt = "  Run extract_stm_geometry.py now? [y/N] "
        try:
            try:
                reply = raw_input(prompt)  # noqa: F821  (Python 2)
            except NameError:
                reply = input(prompt)
        except EOFError:
            print("\n  no input available (stdin is not interactive);")
            print("  stopping. Pass --yes to run the extractor unattended.")
            return 1 if bad else 0
        except KeyboardInterrupt:
            print("\n  interrupted.")
            return 130
        go = reply.strip().lower() in ("y", "yes")

    if not go:
        print("\n  not running the extractor. To run it later:")
        print('    "%s" %s' % (sys.executable, os.path.relpath(EXTRACTOR, ROOT)))
        return 1 if bad else 0

    print("\n" + "=" * 62)
    print("extract_stm_geometry.py")
    print("=" * 62)
    code = subprocess.call([sys.executable, EXTRACTOR], cwd=ROOT)
    if code != 0:
        print("\n  extractor FAILED (exit %d)" % code)
        return code
    print("\n  extractor ok -- CSVs are in %s" % OUTDIR)
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
