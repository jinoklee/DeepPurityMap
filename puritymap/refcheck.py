# =============================================================================
# refcheck.py   Check that a BAM and the reference FASTA are the same genome
#
#   Catches the mismatch before GATK starts, not hours later.
#   Reads without samtools (BAM is BGZF = concatenated gzip, so the gzip module can open it).
# =============================================================================

import gzip
import struct
from pathlib import Path
from typing import Dict, List, Optional, Tuple

# Compared contigs: main chromosomes only (alt · decoy · HLA may differ between BAMs without affecting analysis)
MAIN = [str(i) for i in range(1, 23)] + ["X", "Y"]


class RefError(ValueError):
    pass


def _bam_header_text(bam: Path) -> str:
    """The SAM text header of a BAM (reads the header only)."""
    try:
        with gzip.open(bam, "rb") as f:
            if f.read(4) != b"BAM":
                raise RefError("not a BAM file (magic mismatch)")
            (l_text,) = struct.unpack("<i", f.read(4))
            return f.read(l_text).decode("utf-8", errors="replace")
    except (OSError, EOFError, struct.error) as e:
        raise RefError(f"cannot read the BAM header: {e}")


def bam_read_groups(bam: Path) -> List[str]:
    """Sample names (SM) of the @RG lines in the BAM header. GATK HaplotypeCaller needs at least one @RG with SM."""
    sms = []
    for line in _bam_header_text(bam).splitlines():
        if line.startswith("@RG"):
            sm = [t[3:] for t in line.split("	") if t.startswith("SM:")]
            sms.append(sm[0] if sm else "")
    return sms


def bam_sort_order(bam: Path) -> str:
    """SO of the @HD line ("coordinate", "queryname", "unsorted", "unknown"; "" when absent).
    samtools index · GATK need coordinate-sorted BAMs."""
    for line in _bam_header_text(bam).splitlines():
        if line.startswith("@HD"):
            so = [t[3:] for t in line.split("	") if t.startswith("SO:")]
            return so[0] if so else ""
    return ""


def bam_contigs(bam: Path) -> Dict[str, int]:
    """Contig name → length from the BAM header (reads the header only)."""
    if bam.suffix.lower() == ".cram":
        raise RefError("CRAM is not supported → convert to BAM (samtools view -b)")
    try:
        with gzip.open(bam, "rb") as f:
            if f.read(4) != b"BAM\x01":
                raise RefError("not a BAM file (magic mismatch)")
            (l_text,) = struct.unpack("<i", f.read(4))
            f.read(l_text)                                   # skip the SAM text header
            (n_ref,) = struct.unpack("<i", f.read(4))
            out = {}
            for _ in range(n_ref):
                (l_name,) = struct.unpack("<i", f.read(4))
                name = f.read(l_name)[:-1].decode()           # NUL-terminated
                (l_ref,) = struct.unpack("<i", f.read(4))
                out[name] = l_ref
            return out
    except (OSError, EOFError, struct.error) as e:
        raise RefError(f"cannot read the BAM header: {e}")


def dict_contigs(path: Path) -> Dict[str, int]:
    """Picard .dict (@SQ SN: LN:) or .fai (name\\tlength) → contig lengths."""
    out = {}
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if path.suffix == ".fai":
            f = line.split("\t")
            if len(f) >= 2:
                out[f[0]] = int(f[1])
        elif line.startswith("@SQ"):
            tags = dict(t.split(":", 1) for t in line.split("\t")[1:] if ":" in t)
            if "SN" in tags and "LN" in tags:
                out[tags["SN"]] = int(tags["LN"])
    return out


def bed_problems(bed: Path, ref: Optional[Dict[str, int]] = None, max_lines: int = 5000) -> List[str]:
    """Format problems in a BED (first max_lines data lines). Empty means it looks fine.

    BED = tab-separated, first three columns chrom start end (0-based start < end).
    track/browser/# header lines are allowed. If ref contigs are given, chromosome names must exist there."""
    probs, chroms, n = [], set(), 0
    with open(bed, encoding="utf-8", errors="replace") as f:
        for i, line in enumerate(f, 1):
            line = line.rstrip("\r\n")
            if not line or line.startswith(("track", "browser", "#")):
                continue
            n += 1
            c = line.split("\t")
            if len(c) < 3:
                return [f"line {i}: fewer than 3 tab-separated columns → {line[:60]!r}"]
            if not (c[1].isdigit() and c[2].isdigit()):
                return [f"line {i}: columns 2-3 are not integer start/end (first three columns must be chrom start end) → {line[:60]!r}"]
            if int(c[2]) <= int(c[1]):
                return [f"line {i}: end ≤ start → {line[:60]!r}"]
            chroms.add(c[0])
            if n >= max_lines:
                break
    if n == 0:
        return ["no data lines"]
    if ref:
        missing = sorted(chroms - set(ref))
        if missing:
            probs.append(f"chromosome names not in the reference: {', '.join(missing[:5])}"
                         + (" … (chr prefix?)" if any(not m.lower().startswith("chr") for m in missing) else ""))
    return probs


def _main(contigs: Dict[str, int]) -> Dict[str, Tuple[str, int]]:
    """Main chromosomes only, keyed independently of the 'chr' prefix (1..22, X, Y)."""
    out = {}
    for name, ln in contigs.items():
        k = name[3:] if name.lower().startswith("chr") else name
        if k in MAIN:
            out[k] = (name, ln)
    return out


def compare(bam: Dict[str, int], ref: Dict[str, int]) -> List[str]:
    """List of mismatch descriptions. Empty means they match."""
    b, r = _main(bam), _main(ref)
    if not b:
        return ["BAM header has no main chromosomes (1–22, X, Y)"]
    if not r:
        return ["reference has no main chromosomes (1–22, X, Y)"]
    probs = []
    bn, rn = next(iter(b.values()))[0], next(iter(r.values()))[0]
    if bn.lower().startswith("chr") != rn.lower().startswith("chr"):
        probs.append(f"chromosome naming differs: BAM '{bn}' vs reference '{rn}' (chr prefix)")
    diff = [k for k in b if k in r and b[k][1] != r[k][1]]
    if diff:
        k = diff[0]
        probs.append(f"{len(diff)} chromosome lengths differ (e.g. {k}: BAM {b[k][1]:,} vs reference {r[k][1]:,}) "
                     "→ different genome build (hg19/hg38?)")
    lack = [k for k in b if k not in r]
    if lack:
        probs.append(f"chromosomes missing from the reference: {', '.join(lack)}")
    return probs
