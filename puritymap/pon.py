# =============================================================================
# pon.py   Platform preparation — build the kNN PoN (python -m puritymap prepare-pon)
#
#   Needed only for tumor-only analysis. Once per platform (capture kit).
#   capture BED → tiles + GC (make_tiles.R) → mosdepth --by tiles per normal → build_pon.R
#   Output: PON_DIR/<platform>/  tiles.bed · tiles.gc.tsv · depth/ · normals.tsv · knn_reference.rds · knn_sample_qc.tsv · pon.json
#
#   The sheet has the same format as for analysis (sample_id patient_id type bam platform cancer_type).
#   Only rows with type=normal and a matching platform are used.
# =============================================================================

from pathlib import Path

from . import check, settings as settings_mod, sheet as sheet_mod
from .runner import MAKE_TILES_R, TRACK1, StepError, indexed_bam, newer, run_logged, say

BUILD_PON_R = TRACK1 / "prepare" / "build_pon.R"


def run(conf: str, sheet: str, platform: str, dry: bool = False, rebuild: bool = False) -> int:
    plat = next((p for p in settings_mod.PLATFORMS if p.lower() == platform.lower()), None)
    if plat is None:
        say(f"--platform must be one of {' / '.join(settings_mod.PLATFORMS)}: {platform}")
        return 2

    rep = check.run(conf, sheet, tools=not dry, pon=False)
    errs = [i for i in rep.items if i.level == check.ERR]
    if errs:
        say("pre-run check found errors — fix them first (see everything with python -m puritymap check)")
        for i in errs:
            say(f"  [{i.section}] {i.message}")
        return 1

    st = settings_mod.load(conf)
    if st.pon_dir is None:
        say("PON_DIR is not set → set the folder where the PoN is stored")
        return 1
    sh = sheet_mod.load(sheet, st.bam_dir)
    normals = [s for s in sh.samples if s.type == "normal" and s.platform == plat]
    skipped = [s for s in sh.samples if not (s.type == "normal" and s.platform == plat)]
    if skipped:
        say(f"{len(skipped)} sheet rows skipped (not normal, or platform is not {plat})")
    if not normals:
        say(f"no {plat} normals in the sheet")
        return 1

    out = st.pon_for(plat)
    if rebuild and not dry:
        for f in ("knn_reference.rds", "pon.json", "knn_sample_qc.tsv"):
            (out / f).unlink(missing_ok=True)
    say(f"building PoN — {plat} · {len(normals)} normals ({len({s.patient_id for s in normals})} patients)  → {out}")
    log = out / "logs" / "prepare_pon.log"
    try:
        # 1 tiles + GC
        tiles, tiles_gc = out / "tiles.bed", out / "tiles.gc.tsv"
        if tiles.is_file() and tiles_gc.is_file():
            say("[tiles] exist → reused")
        else:
            say(f"[tiles] {st.bed_for(plat).name} → tiles + GC")
            run_logged([st.tools["RSCRIPT"], MAKE_TILES_R, "--bed", st.bed_for(plat), "--platform", plat, "--out-bed", tiles,
                        "--out-gc", tiles_gc, "--ref-fasta", st.ref_fasta, "--bedtools", st.tools["BEDTOOLS"]], log, dry=dry)

        # 2 normal depth (on tiles — the PoN must use exactly the tile regions)
        depth_dir = out / "depth"
        rows = []
        for i, s in enumerate(normals, 1):
            dep = depth_dir / f"{s.sample_id}.regions.bed.gz"
            rows.append(f"{s.sample_id}\t{s.patient_id}\tdepth/{dep.name}")
            if newer(Path(str(dep) + ".csi"), tiles):
                continue
            say(f"[depth {i}/{len(normals)}] {s.sample_id}")
            bam = indexed_bam(s.bam, s.sample_id, st, log, dry)
            if not dry:
                depth_dir.mkdir(parents=True, exist_ok=True)
            run_logged(["mosdepth", "--by", tiles, "--threads", str(st.threads), "-n",
                        depth_dir / s.sample_id, bam], log, dry=dry)
        table = out / "normals.tsv"
        if not dry:
            table.write_text("sample_id\tpatient_id\tdepth\n" + "\n".join(rows) + "\n", encoding="utf-8")

        # 3 reference
        say("[build] QC → region filter → normalization → GC correction → LOO kNN")
        run_logged([st.tools["RSCRIPT"], BUILD_PON_R, "--normals", table, "--tiles", tiles, "--tiles-gc", tiles_gc,
                    "--platform", plat, "--out", out, "--ref-fasta", st.ref_fasta], log, dry=dry)
    except StepError as e:
        say(f"failed: {e}")
        return 1
    say(f"done → {out}  (per-normal status: knn_sample_qc.tsv)")
    return 0
