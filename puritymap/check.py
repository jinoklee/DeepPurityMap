# =============================================================================
# check.py   Pre-run checks (purity check)
#
#   Checks settings · reference · GATK resources · BED · sheet · BAM · PoN · tools.
#   The goal is to catch problems before GATK runs for hours and then fails.
#   Results are returned as a list of items (the CLI output and the app screen use the same list).
# =============================================================================

import os
import shutil
import subprocess
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import List, Optional

from . import hostpath, refcheck, settings, sheet as sheet_mod

OK, WARN, ERR = "ok", "warn", "error"

# External tools called by Track 1 (all provided inside the Docker image)
TOOLS = ("bash", "gatk", "java", "samtools", "bcftools", "tabix", "bedtools", "mosdepth", "Rscript")


@dataclass
class Item:
    section: str
    level: str
    message: str


class Report:
    def __init__(self):
        self.items: List[Item] = []

    def add(self, section, level, message):
        self.items.append(Item(section, level, hostpath.display(message)))   # show paths as on the user's machine

    @property
    def n_error(self):
        return sum(i.level == ERR for i in self.items)

    @property
    def n_warn(self):
        return sum(i.level == WARN for i in self.items)

    def as_dicts(self):
        return [asdict(i) for i in self.items]


def _file(rep, section, path: Optional[Path], what: str, level_missing=ERR) -> bool:
    if path is None:
        rep.add(section, level_missing, f"{what}: not set")
        return False
    if not path.is_file():
        rep.add(section, level_missing, f"{what} not found: {path}")
        return False
    return True


def _bai(bam: Path) -> Optional[Path]:
    for c in (Path(str(bam) + ".bai"), bam.with_suffix(".bai")):
        if c.is_file():
            return c
    return None


def _writable_dir(path: Path) -> bool:
    p = path
    while not p.exists():                          # if it does not exist yet, the nearest parent it would be created in
        if p.parent == p:
            return False
        p = p.parent
    return p.is_dir() and os.access(p, os.W_OK)


def check_settings(rep: Report, st: settings.Settings):
    S = "settings"
    rep.add(S, OK, f"settings file: {st.path}")
    for k in st.unknown_keys:
        if k == "WGS_BED":
            rep.add(S, WARN, f"WGS_BED cannot be set (WGS always uses {settings.WGS_BED.name}) → ignored")
        else:
            rep.add(S, WARN, f"unknown key '{k}' (typo?) → ignored")
    if _writable_dir(st.results_dir):
        rep.add(S, OK, f"RESULTS_DIR: {st.results_dir}")
    else:
        rep.add(S, ERR, f"cannot write to RESULTS_DIR: {st.results_dir}")
    rep.add(S, OK, f"THREADS={st.threads}  JAVA_MEM={st.java_mem}")

    S = "reference"
    fa = st.ref_fasta
    if _file(rep, S, fa, "REF_FASTA"):
        if fa.suffix == ".gz":
            rep.add(S, ERR, "compressed FASTA (.gz) does not work with bedtools nuc → use the uncompressed file")
        fai = Path(str(fa) + ".fai")
        dct = fa.with_suffix(".dict")
        _file(rep, S, fai, "FASTA index (.fai) — create with samtools faidx")
        _file(rep, S, dct, "FASTA dict (.dict) — create with gatk CreateSequenceDictionary")
        if fai.is_file() and dct.is_file():
            rep.add(S, OK, f"{fa.name} (+ .fai, .dict)")

    S = "GATK resources"
    n_ok = 0
    for key, path in st.resources.items():
        if _file(rep, S, path, f"{key} ({path.name})"):
            if _file(rep, S, Path(str(path) + ".tbi"), f"{key} index ({path.name}.tbi)"):
                n_ok += 1
    if n_ok == len(st.resources):
        rep.add(S, OK, f"all 4 files + .tbi found: {st.gatk_resource_dir}")


def check_sheet(rep: Report, st: settings.Settings, sh: sheet_mod.Sheet, pon: bool = True):
    S = "sample sheet"
    for e in sh.errors:
        rep.add(S, ERR, e)
    for w in sh.warnings:
        rep.add(S, WARN, w)
    jobs = sh.jobs()
    n_pair = sum(j.mode == "paired" for j in jobs)
    rep.add(S, OK, f"{sh.path.name}: {len(sh.samples)} samples → {len(jobs)} analyses "
                   f"(paired {n_pair}, tumor-only {len(jobs) - n_pair})")

    ref_contigs = None
    dct = st.ref_fasta.with_suffix(".dict")
    if dct.is_file():
        ref_contigs = refcheck.dict_contigs(dct)
    elif Path(str(st.ref_fasta) + ".fai").is_file():
        ref_contigs = refcheck.dict_contigs(Path(str(st.ref_fasta) + ".fai"))

    S = "BED"
    for plat in sorted({s.platform for s in sh.samples}):
        bed = st.bed_for(plat)
        if bed is None:
            key = settings.USER_BED_KEYS[plat]
            rep.add(S, ERR, f"{plat} samples present → set {key} (capture kit BED) in the settings file")
        elif _file(rep, S, bed, f"{plat} BED"):
            probs = refcheck.bed_problems(bed, ref_contigs)
            for p in probs:
                rep.add(S, ERR, f"{plat} BED {bed.name}: {p}")
            if not probs:
                rep.add(S, OK, f"{plat}: {bed.name}" + (" (fixed)" if plat == "WGS" else ""))

    S = "BAM"
    for s in sh.samples:
        tag = f"{s.sample_id} ({s.type})"
        if not s.bam.is_file():
            rep.add(S, ERR, f"{tag}: BAM not found → {s.bam}")
            continue
        bai = _bai(s.bam)
        if bai is None:
            rep.add(S, WARN, f"{tag}: no index (.bai) → created automatically at run time")
        elif bai.stat().st_mtime < s.bam.stat().st_mtime:
            rep.add(S, WARN, f"{tag}: .bai is older than the BAM → recreating it is recommended")
        try:
            probs = refcheck.compare(refcheck.bam_contigs(s.bam), ref_contigs) if ref_contigs else None
            sms = refcheck.bam_read_groups(s.bam)
            so = refcheck.bam_sort_order(s.bam)
        except refcheck.RefError as e:
            rep.add(S, ERR, f"{tag}: {e}")
            continue
        if so != "coordinate":                            # samtools index (and so GATK) needs a coordinate-sorted BAM
            out = s.bam.with_name(s.bam.stem + ".sorted.bam").name
            rep.add(S, ERR, f"{tag}: not coordinate-sorted (@HD SO:{so or 'missing'}) → samtools sort -@ 8 -o {out} {s.bam.name} "
                            f"&& samtools index {out}, then use {out} in the sheet; see docs/setup_guide.md")
        if not sms:                                       # HaplotypeCaller stops at once without read groups
            rep.add(S, ERR, f"{tag}: no read group (@RG) in the BAM header → add one with gatk AddOrReplaceReadGroups "
                            f"(RGSM = sample name); see docs/setup_guide.md")
        elif not all(sms):
            rep.add(S, ERR, f"{tag}: @RG without SM (sample name) → set it with gatk AddOrReplaceReadGroups --RGSM")
        if probs is None:
            rep.add(S, WARN, f"{tag}: reference has no index → genome match check skipped")
        elif probs:
            for p in probs:
                rep.add(S, ERR, f"{tag}: {p}")
        else:
            rep.add(S, OK, f"{tag}: {s.bam.name} — matches the reference")

    if not pon:                                    # prepare-pon: the PoN is about to be built, so do not check it
        return
    S = "PoN"
    need = sorted({j.tumor.platform for j in jobs if j.mode == "tumor-only"})
    if not need:
        if jobs:
            rep.add(S, OK, "no tumor-only samples → PoN not needed")
        return
    if st.pon_dir is None:
        rep.add(S, ERR, f"tumor-only samples present ({', '.join(need)}) → set PON_DIR in the settings file and run prepare-pon")
        return
    for plat in need:
        d = st.pon_for(plat)
        ok = all(_file(rep, S, d / f, f"{plat} PoN {f}") for f in ("knn_reference.rds", "pon.json"))
        if ok:
            rep.add(S, OK, f"{plat}: {d}")
        else:
            rep.add(S, ERR, f"{plat}: no PoN → build one from {plat} normals with prepare-pon")


TRACK1_R = ("data.table", "matrixStats", "vcfR", "PSCBS", "DNAcopy", "aroma.light", "MASS", "ggplot2", "viridis", "gtable")


def check_tools(rep: Report, st: settings.Settings):
    S = "tools"
    named = {"gatk": st.tools["GATK"], "Rscript": st.tools["RSCRIPT"], "bedtools": st.tools["BEDTOOLS"]}
    miss = []
    for t in TOOLS:
        exe = shutil.which(named.get(t, t))
        if exe:
            rep.add(S, OK, f"{t}: {exe}")
        else:
            miss.append(t)
    if miss:
        rep.add(S, ERR, f"not found: {', '.join(miss)} → Track 1 cannot run in this environment (install them — see the README; set GATK= etc. in puritymap.conf)")

    # R packages of Track 1 (features · segmentation); PSCBS is installed separately, so it is the one most often missing
    code = "for (p in c(%s)) if (!requireNamespace(p, quietly = TRUE)) cat(p, '')" % ", ".join(f"'{p}'" for p in TRACK1_R)
    try:
        r = subprocess.run([st.tools["RSCRIPT"], "-e", code], capture_output=True, text=True, timeout=120)
        gone = r.stdout.split()
        rep.add(S, ERR if gone or r.returncode else OK,
                f"R packages missing: {' '.join(gone)}" + (" → install.packages(\"PSCBS\")" if "PSCBS" in gone else "") if gone else
                ("R packages: " + " · ".join(TRACK1_R)) if r.returncode == 0 else f"R check failed: {r.stderr.strip()[-200:]}")
    except (OSError, subprocess.TimeoutExpired) as e:
        rep.add(S, ERR, f"cannot run Rscript ({st.tools['RSCRIPT']}): {e}")

    # Prediction environment: check that the packages actually import, not just that the executable exists
    py = st.tools["PREDICT_PYTHON"]
    try:
        r = subprocess.run([py, "-c", "import torch; torch.cuda.is_available() or setattr(torch.cuda, 'get_device_capability', lambda *a, **k: (0, 0)); "
                                      "import natten, transformers; print(torch.__version__, natten.__version__, "
                                      "transformers.__version__, torch.cuda.is_available())"],
                           capture_output=True, text=True, timeout=120)
        if r.returncode == 0:
            tv, nv, hv, cuda = r.stdout.split()
            rep.add(S, OK, f"prediction env: torch {tv} · natten {nv} · transformers {hv} · GPU {'available' if cuda == 'True' else 'none → CPU'}")
        else:
            last = (r.stderr.strip().splitlines() or ["?"])[-1]
            rep.add(S, ERR, f"prediction env ({py}) failed to import packages: {last}")
    except (OSError, subprocess.TimeoutExpired) as e:
        rep.add(S, ERR, f"cannot run prediction env ({py}): {e}")


def run(conf_path, sheet_path=None, tools=True, pon=True) -> Report:
    rep = Report()
    try:
        st = settings.load(conf_path)
    except settings.ConfError as e:
        rep.add("settings", ERR, str(e))
        return rep
    check_settings(rep, st)
    if sheet_path:
        try:
            check_sheet(rep, st, sheet_mod.load(sheet_path, st.bam_dir), pon=pon)
        except sheet_mod.SheetError as e:
            rep.add("sample sheet", ERR, str(e))
    if tools:
        check_tools(rep, st)
    return rep
