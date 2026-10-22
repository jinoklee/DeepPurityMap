# =============================================================================
# runner.py   Run the analysis (python -m puritymap run) for every tumor sample
#
#   1_genotype_depth.sh → 3_features.R → 4_predict.py → 4_purity_figure.R   (genotype_depth · features · purity · figure)
#   - Finished steps are skipped. A step is rerun if an earlier step's output is newer (modification time).
#   - logs/<step>.log per step, metadata.json per sample (relative paths only).
#   - If platform · mode · normal differ from the previous run, that sample stops (so results from different conditions never mix).
#   - If one sample fails, the next sample still runs (change with --stop-on-error).
# =============================================================================

import hashlib
import json
import os
import shutil
import subprocess
import time
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional

from . import __version__, check, hostpath, settings as settings_mod, sheet as sheet_mod
from .settings import ROOT, Settings

TRACK1 = ROOT / "track1"
GENOTYPE_DEPTH_SH = TRACK1 / "1_genotype_depth.sh"
MAKE_TILES_R = TRACK1 / "prepare" / "make_tiles.R"
FEATURES_R = TRACK1 / "3_features.R"
PREDICT_PY = TRACK1 / "4_predict.py"
FIGURE_R = TRACK1 / "4_purity_figure.R"

STEPS = ("genotype_depth", "features", "purity", "figure")


class StepError(RuntimeError):
    pass


def now() -> str:
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")


def say(msg: str):
    print(hostpath.display(msg), flush=True)


def run_logged(cmd: List[str], log: Path, env: Optional[Dict[str, str]] = None, dry: bool = False) -> None:
    """Run a command with its output going to log. On failure raise StepError with the tail of the log.

    The working directory is the folder of log, so files that tools leave in the current directory
    (e.g. lightning_logs/ from GATK NVScoreVariants) end up in that sample's logs/, not in the code folder.
    Commands receive absolute paths only, so changing the working directory has no other effect."""
    shown = " ".join(str(c) for c in cmd)
    if dry:
        say(f"    $ {shown}")
        return
    log.parent.mkdir(parents=True, exist_ok=True)
    with open(log, "a", encoding="utf-8") as f:
        f.write(f"\n# {now()}  $ {shown}\n")
        f.flush()
        r = subprocess.run([str(c) for c in cmd], stdout=f, stderr=subprocess.STDOUT,
                           env={**os.environ, **(env or {})}, cwd=log.parent)
    if r.returncode != 0:
        tail = log.read_text(encoding="utf-8", errors="replace").strip().splitlines()[-15:]
        raise StepError(f"exit code {r.returncode} — log {log}\n      " + "\n      ".join(tail))


def newer(out: Path, *inputs: Path) -> bool:
    """True if out exists and is newer than every input (= no need to redo)."""
    if not out.exists():
        return False
    t = out.stat().st_mtime
    return all(not i.exists() or i.stat().st_mtime <= t for i in inputs)


# ── BAM index: the user's BAM folder is read-only → if missing, link + index in the results folder ──
def indexed_bam(bam: Path, sample_id: str, st: Settings, log: Path, dry: bool) -> Path:
    for bai in (Path(str(bam) + ".bai"), bam.with_suffix(".bai")):
        if bai.is_file() and bai.stat().st_mtime >= bam.stat().st_mtime:
            return bam
    d = st.results_dir / "_bam_index"
    link = d / f"{sample_id}.bam"
    if not dry:
        d.mkdir(parents=True, exist_ok=True)
        if link.is_symlink() or link.exists():
            if link.resolve() != bam.resolve():
                link.unlink()
        if not link.exists():
            link.symlink_to(bam)
    if not newer(Path(str(link) + ".bai"), bam):
        say(f"  · creating BAM index: {bam.name} → {d}")
        run_logged(["samtools", "index", "-@", str(st.threads), link], log, dry=dry)
    return link


# ── tiles (reference regions for paired): built once per capture BED and cached in the results folder ──
def tiles_for(platform: str, st: Settings, dry: bool) -> Dict[str, Path]:
    bed = st.bed_for(platform)
    key = hashlib.md5(bed.read_bytes()).hexdigest()[:8]
    d = st.results_dir / "_tiles" / f"{platform}-{key}"
    out = {"bed": d / "tiles.bed", "gc": d / "tiles.gc.tsv"}
    if not (out["bed"].is_file() and out["gc"].is_file()):
        say(f"[tiles] {platform}: {bed.name} → tiles + GC  ({d})")
        run_logged([st.tools["RSCRIPT"], MAKE_TILES_R, "--bed", bed, "--platform", platform, "--out-bed", out["bed"], "--out-gc", out["gc"],
                    "--ref-fasta", st.ref_fasta, "--bedtools", st.tools["BEDTOOLS"]],
                   st.results_dir / "_tiles" / f"{platform}-{key}.log", dry=dry)
    return out


# ── metadata.json ────────────────────────────────────────────────────────────
def load_meta(path: Path) -> dict:
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return {}


def save_meta(path: Path, meta: dict):
    tmp = path.with_suffix(".json.tmp")
    tmp.write_text(json.dumps(meta, ensure_ascii=False, indent=2), encoding="utf-8")   # utf-8 required (avoid Windows cp949)
    tmp.replace(path)


def run_job(job: sheet_mod.Job, st: Settings, redo: set, dry: bool) -> dict:
    t, n = job.tumor, job.normal
    sid = t.sample_id
    sdir = st.results_dir / sid
    logs = sdir / "logs"
    meta_path = sdir / "metadata.json"
    meta = load_meta(meta_path)

    ident = {"platform": t.platform, "mode": job.mode, "normal_id": n.sample_id if n else None}
    old = {k: meta.get(k) for k in ident} if meta else None
    if old and any(old[k] != ident[k] for k in ident) and "all" not in redo:
        raise StepError(f"conditions differ from the previous run (before {old} / now {ident}) → use a different sample_id or --redo all")
    if "all" in redo and not dry and sdir.exists():
        for sub in ("1_genotype", "2_depth", "3_features", "4_purity"):
            shutil.rmtree(sdir / sub, ignore_errors=True)
        meta = {}

    meta.update({"sample_id": sid, "patient_id": t.patient_id, "cancer_type": t.cancer_type, **ident,
                 "tumor_bam": hostpath.display(str(t.bam)), "normal_bam": hostpath.display(str(n.bam)) if n else None,
                 "reference": st.ref_fasta.name, "pipeline_version": __version__})
    meta.setdefault("steps", {})

    def step(name: str, done: bool, fn):
        rec = meta["steps"].get(name, {})
        if done and name not in redo:
            say(f"  [{name}] skipped (already done)")
            return
        say(f"  [{name}] start")
        t0 = time.time()
        rec.update({"status": "running", "started": now(), "log": f"logs/{name}.log"})
        meta["steps"][name] = rec
        if not dry:
            sdir.mkdir(parents=True, exist_ok=True)
            save_meta(meta_path, meta)
        try:
            fn(logs / f"{name}.log")
        except StepError:
            rec.update({"status": "failed", "ended": now(), "elapsed_s": round(time.time() - t0)})
            if not dry:
                save_meta(meta_path, meta)
            raise
        rec.update({"status": "success" if not dry else "dry-run", "ended": now(), "elapsed_s": round(time.time() - t0)})
        if not dry:
            save_meta(meta_path, meta)
        say(f"  [{name}] done ({rec['elapsed_s']}s)")

    # ── 1·2 genotype · depth ────────────────────────────────────────────────
    g_vcf = sdir / "1_genotype" / "final.vcf.gz"
    t_dep = sdir / "2_depth" / f"{sid}.regions.bed.gz"
    n_dep = sdir / "2_depth" / f"{n.sample_id}.regions.bed.gz" if n else None
    gd_done = Path(str(g_vcf) + ".tbi").exists() and Path(str(t_dep) + ".csi").exists() and \
        (n_dep is None or Path(str(n_dep) + ".csi").exists())

    def do_gd(log):
        if "genotype_depth" in redo and not dry:
            for sub in ("1_genotype", "2_depth"):
                shutil.rmtree(sdir / sub, ignore_errors=True)
        tb = indexed_bam(t.bam, sid, st, log, dry)
        cmd = ["bash", GENOTYPE_DEPTH_SH, "-t", tb, "-s", sid, "-b", st.bed_for(t.platform), "-o", sdir]
        if n:
            cmd += ["-n", indexed_bam(n.bam, n.sample_id, st, log, dry), "-N", n.sample_id]
        if t.platform == "WGS":                            # WGS: depth on the fixed 1kb windows for tumor and normal
            cmd += ["-d", settings_mod.WGS_BED]
        run_logged(cmd, log, env=st.as_env(), dry=dry)
    step("genotype_depth", gd_done, do_gd)

    # ── 3 features ──────────────────────────────────────────────────────────
    fdir = sdir / "3_features"
    seg_png = fdir / f"{sid}.seg.png"
    f_done = newer(seg_png, g_vcf, t_dep, *([n_dep] if n_dep else [])) and (fdir / "qc.tsv").exists()

    def do_features(log):
        cmd = [st.tools["RSCRIPT"], FEATURES_R, "--sample", sid, "--platform", t.platform, "--vcf", g_vcf, "--tumor-depth", t_dep, "--out", fdir]
        if n:
            tl = tiles_for(t.platform, st, dry)
            cmd += ["--normal-depth", n_dep, "--tiles", tl["bed"], "--tiles-gc", tl["gc"]]
            meta["tiles"] = Path(os.path.relpath(tl["bed"], sdir)).as_posix()     # relative path inside the results folder
        else:
            pon = st.pon_for(t.platform)
            pj = load_meta(pon / "pon.json")
            if pj.get("reference") and pj["reference"] != st.ref_fasta.name:
                raise StepError(f"PoN reference ({pj['reference']}) differs from the current reference ({st.ref_fasta.name})")
            cmd += ["--pon", pon / "knn_reference.rds", "--patient", t.patient_id]
            meta["pon"] = {"dir": hostpath.display(str(pon)), "created": pj.get("created"),
                           "n_normals": pj.get("n_normals_used"), "tiles_md5": pj.get("tiles_md5")}
        run_logged(cmd, log, dry=dry)
    step("features", f_done, do_features)

    # ── 4 purity ────────────────────────────────────────────────────────────
    pdir = sdir / "4_purity"
    p_tsv = pdir / f"{sid}.purity.tsv"

    def do_purity(log):
        run_logged([st.tools["PREDICT_PYTHON"], PREDICT_PY, "--sample", sid, "--features-dir", fdir, "--out", pdir],
                   log, dry=dry)
    step("purity", newer(p_tsv, seg_png), do_purity)

    if p_tsv.exists():
        head, row = (l.split("\t") for l in p_tsv.read_text(encoding="utf-8").splitlines()[:2])
        rec = dict(zip(head, row))                         # by column name, so added columns never break this
        meta["purity"] = {k: rec[k] for k in ("model", "model_version", "device", "weights_sha256") if k in rec}
        meta["purity"]["predicted_purity"] = float(rec["predicted_purity"])
        if not dry:
            save_meta(meta_path, meta)

    # ── purity figure: SNP depth–VAF distribution with the expected positions for each purity ──
    seg_txt = fdir / f"{sid}.pscbs_seg.txt"
    feat_txt = fdir / f"{sid}.feature.txt"
    fig_png = pdir / f"{sid}_purity.png"

    def do_figure(log):
        run_logged([st.tools["RSCRIPT"], FIGURE_R, "--sample", sid, "--seg", seg_txt, "--features", feat_txt,
                    "--out", pdir], log, dry=dry)
    step("figure", newer(fig_png, seg_txt, feat_txt), do_figure)
    return meta


def run(conf: str, sheet: str, only: Optional[List[str]] = None, redo: Optional[List[str]] = None,
        dry: bool = False, stop_on_error: bool = False) -> int:
    redo = set(redo or [])
    bad = redo - set(STEPS) - {"all"}
    if bad:
        say(f"--redo must be from {', '.join(STEPS)}, all: {', '.join(bad)}")
        return 2

    rep = check.run(conf, sheet, tools=not dry)
    errs = [i for i in rep.items if i.level == check.ERR]
    if errs:
        say("pre-run check found errors — fix them first (see everything with python -m puritymap check)")
        for i in errs:
            say(f"  [{i.section}] {i.message}")
        return 1

    st = settings_mod.load(conf)
    jobs = sheet_mod.load(sheet, st.bam_dir).jobs()
    if only:
        miss = set(only) - {j.tumor.sample_id for j in jobs}
        if miss:
            say(f"--samples has tumors not in the sheet: {', '.join(sorted(miss))}")
            return 2
        jobs = [j for j in jobs if j.tumor.sample_id in only]
    say(f"Analysis — {len(jobs)} samples{' (dry-run: commands only)' if dry else ''}  → {st.results_dir}")

    results = []
    for i, job in enumerate(jobs, 1):
        say(f"\n== [{i}/{len(jobs)}] {job.tumor.sample_id}  ({job.tumor.platform}, {job.mode})")
        try:
            meta = run_job(job, st, redo, dry)
            pur = meta.get("purity", {}).get("predicted_purity")
            results.append((job, "ok" if not dry else "dry-run", pur))
        except (StepError, OSError) as e:
            say(f"  failed: {e}")
            results.append((job, "failed", None))
            if stop_on_error:
                break

    say("\n" + "─" * 60)
    for job, status, pur in results:
        say(f"  {job.tumor.sample_id:<30} {job.mode:<11} {status:<6} {'' if pur is None else f'purity {pur:.4f}'}")
    n_fail = sum(s == "failed" for _, s, _ in results)
    say(f"done {len(results) - n_fail} · failed {n_fail}")
    return 1 if n_fail else 0
