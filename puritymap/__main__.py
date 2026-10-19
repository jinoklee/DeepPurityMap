# =============================================================================
# python -m puritymap <command>
#
#   check        pre-run check    python -m puritymap check --conf puritymap.conf --sheet samples.tsv
#   run          run the analysis python -m puritymap run --conf puritymap.conf --sheet samples.tsv
#   prepare-pon  build the PoN    python -m puritymap prepare-pon --conf puritymap.conf --sheet normals.tsv --platform WES
#
#   --conf · --sheet take the user's machine paths as is (mapped automatically inside the container).
# =============================================================================

import argparse
import json
import sys
from pathlib import Path

from . import __version__, check, pon, runner
from .settings import resolve_path

MARK = {check.OK: "  OK   ", check.WARN: "  WARN ", check.ERR: "  ERROR"}


def cmd_check(args) -> int:
    rep = check.run(args.conf, args.sheet, tools=not args.skip_tools)
    if args.json:
        print(json.dumps({"errors": rep.n_error, "warnings": rep.n_warn, "items": rep.as_dicts()},
                         ensure_ascii=False, indent=2))
        return 1 if rep.n_error else 0

    section = None
    for it in rep.items:
        if it.section != section:
            section = it.section
            print(f"\n[{section}]")
        print(f"{MARK[it.level]}  {it.message}")
    print(f"\nresult: {rep.n_error} errors · {rep.n_warn} warnings"
          + ("  → fix the errors and check again" if rep.n_error else "  → ready to run"))
    return 1 if rep.n_error else 0


def cmd_run(args) -> int:
    only = [s.strip() for s in args.samples.split(",")] if args.samples else None
    redo = [s.strip() for s in args.redo.split(",")] if args.redo else None
    return runner.run(args.conf, args.sheet, only=only, redo=redo, dry=args.dry_run,
                      stop_on_error=args.stop_on_error)


def cmd_pon(args) -> int:
    return pon.run(args.conf, args.sheet, args.platform, dry=args.dry_run, rebuild=args.rebuild)


def main(argv=None) -> int:
    # Keep the Windows console (cp949) from garbling or crashing on non-ASCII symbols
    for s in (sys.stdout, sys.stderr):
        try:
            s.reconfigure(encoding="utf-8", errors="replace")
        except (AttributeError, ValueError):
            pass

    ap = argparse.ArgumentParser(prog="puritymap", description=f"DeepPurityMap pipeline {__version__}")
    ap.add_argument("--version", action="version", version=__version__)
    sub = ap.add_subparsers(dest="cmd", required=True)

    def common(p, sheet_required=True):
        p.add_argument("--conf", required=True, help="settings file (puritymap.conf)")
        p.add_argument("--sheet", required=sheet_required, help="sample sheet (tab-separated)")

    p = sub.add_parser("check", help="pre-run check (settings · reference · sheet · BAM · PoN · tools)")
    common(p, sheet_required=False)
    p.add_argument("--skip-tools", action="store_true", help="skip the external tool check")
    p.add_argument("--json", action="store_true", help="output JSON")
    p.set_defaults(func=cmd_check)

    p = sub.add_parser("run", help="run the analysis: genotype · depth → features → purity → purity figure")
    common(p)
    p.add_argument("--samples", help="subset: tumor sample_ids, comma-separated")
    p.add_argument("--redo", help="steps to redo: genotype_depth, features, purity, figure, all (comma-separated)")
    p.add_argument("--dry-run", action="store_true", help="show the commands without running them")
    p.add_argument("--stop-on-error", action="store_true", help="stop when a sample fails")
    p.set_defaults(func=cmd_run)

    p = sub.add_parser("prepare-pon", help="platform preparation: build a kNN PoN from normals (for tumor-only)")
    common(p)
    p.add_argument("--platform", required=True, help="WES / Panel / WGS")
    p.add_argument("--rebuild", action="store_true", help="delete the existing PoN and rebuild (tiles · depth are reused)")
    p.add_argument("--dry-run", action="store_true", help="show the commands without running them")
    p.set_defaults(func=cmd_pon)

    args = ap.parse_args(argv)
    if getattr(args, "conf", None):
        args.conf = str(resolve_path(args.conf, Path.cwd()))
    if getattr(args, "sheet", None):
        args.sheet = str(resolve_path(args.sheet, Path.cwd()))
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
