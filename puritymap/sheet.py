# =============================================================================
# sheet.py   Read the sample sheet (samples.tsv)
#
#   Columns: sample_id  patient_id  type  bam  platform  cancer_type
#   type is given explicitly by the user (never guessed from file names).
#   Mode is automatic: paired if a normal with the same patient_id exists, otherwise tumor-only.
#   Relative bam paths are resolved against BAM_DIR of the settings file when set, else the folder that holds the sheet.
# =============================================================================

import csv
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional

from .settings import PLATFORMS, resolve_path

COLUMNS = ("sample_id", "patient_id", "type", "bam", "platform", "cancer_type")
TYPES = ("tumor", "normal")
SAMPLE_ID_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_-]*$")    # restricted because it becomes a result folder name


class SheetError(ValueError):
    pass


@dataclass
class Sample:
    line: int                    # line number in the sheet (for messages)
    sample_id: str
    patient_id: str
    type: str
    bam: Path
    platform: str
    cancer_type: str


@dataclass
class Job:
    """One analysis unit per tumor. Paired if a normal exists."""
    tumor: Sample
    normal: Optional[Sample]

    @property
    def mode(self) -> str:
        return "paired" if self.normal else "tumor-only"


@dataclass
class Sheet:
    path: Path
    samples: List[Sample]
    errors: List[str]
    warnings: List[str]

    def jobs(self) -> List[Job]:
        normals: Dict[str, Sample] = {}
        for s in self.samples:
            if s.type == "normal":
                normals.setdefault(s.patient_id, s)
        return [Job(s, normals.get(s.patient_id)) for s in self.samples if s.type == "tumor"]


def _norm_platform(v: str) -> Optional[str]:
    for p in PLATFORMS:                               # case-insensitive (wes, PANEL …)
        if v.lower() == p.lower():
            return p
    return None


def load(sheet_path, bam_dir=None) -> Sheet:
    sheet_path = Path(sheet_path).resolve()
    if not sheet_path.is_file():
        raise SheetError(f"sample sheet not found: {sheet_path}")
    text = sheet_path.read_text(encoding="utf-8-sig")       # strip the BOM that Excel adds
    rows = list(csv.reader(text.splitlines(), delimiter="\t"))
    rows = [(i, r) for i, r in enumerate(rows, 1) if any(c.strip() for c in r) and not r[0].startswith("#")]
    if not rows:
        raise SheetError(f"sample sheet is empty: {sheet_path}")

    hl, header = rows[0]
    header = [h.strip().lower() for h in header]
    if len(header) == 1 and "," in header[0]:
        raise SheetError("columns are not tab-separated (commas?) → in Excel save as 'Text (Tab delimited)'")
    miss = [c for c in COLUMNS if c not in header]
    if miss:
        raise SheetError(f"line {hl}: header is missing columns: {', '.join(miss)}  (required: {' '.join(COLUMNS)})")
    idx = {c: header.index(c) for c in COLUMNS}

    errors, warnings, samples = [], [], []
    base = Path(bam_dir) if bam_dir else sheet_path.parent
    for n, r in rows[1:]:
        get = lambda c: r[idx[c]].strip() if idx[c] < len(r) else ""
        sid, pid, typ = get("sample_id"), get("patient_id"), get("type").lower()
        bam, plat, cancer = get("bam"), _norm_platform(get("platform")), get("cancer_type")
        where = f"line {n} ({sid or '?'})"

        bad = False
        if not SAMPLE_ID_RE.match(sid):
            errors.append(f"{where}: sample_id may contain only letters, digits, '-' and '_' (must start with a letter or digit)")
            bad = True
        if not pid:
            errors.append(f"{where}: patient_id is empty")
            bad = True
        if typ not in TYPES:
            errors.append(f"{where}: type must be tumor or normal → {get('type')!r}")
            bad = True
        if not bam:
            errors.append(f"{where}: bam path is empty")
            bad = True
        if plat is None:
            errors.append(f"{where}: platform must be one of {' / '.join(PLATFORMS)} → {get('platform')!r}")
            bad = True
        if bad:
            continue
        if typ == "tumor" and not cancer:
            warnings.append(f"{where}: cancer_type is empty → shown as 'not specified' in the report")
        try:
            bam_path = resolve_path(bam, base)
        except ValueError as e:
            errors.append(f"{where}: {e}")
            continue
        samples.append(Sample(n, sid, pid, typ, bam_path, plat, cancer))

    # ── Sheet-wide rules ──────────────────────────────────────────────────────
    seen: Dict[str, int] = {}
    for s in samples:
        if s.sample_id in seen:
            errors.append(f"line {s.line}: duplicate sample_id '{s.sample_id}' (same as line {seen[s.sample_id]})")
        seen.setdefault(s.sample_id, s.line)

    bams: Dict[Path, str] = {}
    for s in samples:
        if s.bam in bams:
            errors.append(f"line {s.line}: same BAM is also used by '{bams[s.bam]}' → {s.bam}")
        bams.setdefault(s.bam, s.sample_id)

    by_patient: Dict[str, List[Sample]] = {}
    for s in samples:
        by_patient.setdefault(s.patient_id, []).append(s)
    for pid, ss in by_patient.items():
        nor = [s for s in ss if s.type == "normal"]
        tum = [s for s in ss if s.type == "tumor"]
        if len(nor) > 1:
            errors.append(f"patient_id '{pid}': {len(nor)} normals ({', '.join(s.sample_id for s in nor)}) "
                          "→ keep only one matched normal")
        if nor and not tum:
            warnings.append(f"patient_id '{pid}': normal only, no tumor → not analyzed "
                            "(for a PoN, put it in the prepare-pon sheet)")
        for t in tum:
            for nm in nor:
                if nm.platform != t.platform:
                    errors.append(f"patient_id '{pid}': tumor ({t.platform}) and normal ({nm.platform}) platforms differ")

    if not any(s.type == "tumor" for s in samples) and not errors:
        warnings.append("no tumor rows → nothing to analyze")

    return Sheet(sheet_path, samples, errors, warnings)
