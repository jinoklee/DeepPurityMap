# =============================================================================
# settings.py   Read the settings file (puritymap.conf)
#
#   Format: one KEY=VALUE per line. Anything after '#' is a comment. Quotes around values are stripped.
#   Relative paths are resolved against the folder that holds the settings file (not the working directory).
#   The format can be sourced by bash and used with docker compose --env-file.
# =============================================================================

import os
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional

from . import hostpath

# Repository root (parent of puritymap/). Bundled files such as the default BED are found from here.
ROOT = Path(__file__).resolve().parent.parent

PLATFORMS = ("WES", "Panel", "WGS")

# GATK resources — Broad standard file names inside GATK_RESOURCE_DIR. Each can be overridden by its key.
GATK_RESOURCES = {
    "DBSNP":        "dbsnp_146.hg38.vcf.gz",
    "MILLS":        "Mills_and_1000G_gold_standard.indels.hg38.vcf.gz",
    "KNOWN_INDELS": "Homo_sapiens_assembly38.known_indels.vcf.gz",
    "OMNI":         "1000G_omni2.5.hg38.vcf.gz",
}

# GATK intervals / capture BED per platform
#   WGS        fixed (bundled). Users cannot change it.
#   WES/Panel  differ by capture kit, so the user sets them in the settings file (no default).
WGS_BED = ROOT / "resources" / "bed" / "WGS_snp.1kb.bed"
USER_BED_KEYS = {"WES": "WES_BED", "Panel": "PANEL_BED"}

# Tool locations. Inside the image they are preset via environment variables (PURITY_*), so no setting is needed.
# Set them in the settings file only when running directly on a server without Docker (command name or path).
TOOL_KEYS = {
    "GATK":           ("PURITYMAP_GATK", "gatk"),
    "RSCRIPT":        ("PURITYMAP_RSCRIPT", "Rscript"),
    "PREDICT_PYTHON": ("PURITYMAP_PREDICT_PYTHON", sys.executable or "python3"),
    "BEDTOOLS":       ("PURITYMAP_BEDTOOLS", "bedtools"),
}

KNOWN_KEYS = {"REF_FASTA", "GATK_RESOURCE_DIR", "RESULTS_DIR", "PON_DIR", "BAM_DIR",
              "THREADS", "JAVA_MEM", "NCBI_API_KEY", *USER_BED_KEYS.values(), *GATK_RESOURCES, *TOOL_KEYS}

REQUIRED_KEYS = ("REF_FASTA", "GATK_RESOURCE_DIR", "RESULTS_DIR")


class ConfError(ValueError):
    pass


@dataclass
class Settings:
    path: Path                                    # the settings file itself
    ref_fasta: Path
    gatk_resource_dir: Path
    results_dir: Path
    pon_dir: Optional[Path]
    threads: int
    java_mem: str
    beds: Dict[str, Optional[Path]]               # platform → intervals BED (None if not set)
    resources: Dict[str, Path]                    # DBSNP … → file
    tools: Dict[str, str] = field(default_factory=dict)   # GATK · RSCRIPT · PREDICT_PYTHON · BEDTOOLS
    bam_dir: Optional[Path] = None                # relative bam paths in the sheet are read from here
    ncbi_api_key: Optional[str] = None            # PubMed (report); personal key, from the settings file or the environment
    unknown_keys: List[str] = field(default_factory=list)

    def bed_for(self, platform: str) -> Optional[Path]:
        return self.beds.get(platform)

    def pon_for(self, platform: str) -> Optional[Path]:
        return self.pon_dir / platform if self.pon_dir else None

    def as_env(self) -> Dict[str, str]:
        """Environment variables passed to subprocesses (sh · R). All paths are absolute."""
        env = {
            "REF_FASTA": str(self.ref_fasta),
            "RESULTS_DIR": str(self.results_dir),
            "THREADS": str(self.threads),
            "JAVA_MEM": self.java_mem,
            "GATK": self.tools["GATK"],
            **{k: str(v) for k, v in self.resources.items()},
        }
        if self.pon_dir:
            env["PON_DIR"] = str(self.pon_dir)
        if self.ncbi_api_key:
            env["NCBI_API_KEY"] = self.ncbi_api_key
        return env


def parse_conf(text: str) -> Dict[str, str]:
    out = {}
    for n, raw in enumerate(text.splitlines(), 1):
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        if line.startswith("export "):
            line = line[len("export "):].lstrip()
        if "=" not in line:
            raise ConfError(f"line {n}: not in KEY=VALUE form → {raw!r}")
        key, val = line.split("=", 1)
        key, val = key.strip(), val.strip()
        if val[:1] in ("'", '"'):                     # quoted value: up to the closing quote
            q = val[0]
            end = val.find(q, 1)
            if end < 0:
                raise ConfError(f"line {n}: unclosed quote → {raw!r}")
            val = val[1:end]
        else:                                          # unquoted value: ' #' starts a comment
            i = val.find(" #")
            if i >= 0:
                val = val[:i].rstrip()
        if not key.replace("_", "").isalnum():
            raise ConfError(f"line {n}: invalid key name → {key!r}")
        out[key] = val
    return out


def resolve_path(val: str, base: Path) -> Path:
    """Path as written by the user → openable absolute path. Inside the container it maps under /host; relative paths use base."""
    p = Path(hostpath.to_local(val.strip()).replace("\\", "/")).expanduser()
    return (p if p.is_absolute() else base / p).resolve()


def load(conf_path) -> Settings:
    conf_path = Path(conf_path).resolve()
    if not conf_path.is_file():
        raise ConfError(f"settings file not found: {conf_path}")
    kv = parse_conf(conf_path.read_text(encoding="utf-8-sig"))
    base = conf_path.parent

    missing = [k for k in REQUIRED_KEYS if not kv.get(k)]
    if missing:
        raise ConfError(f"required keys missing: {', '.join(missing)}  ({conf_path.name})")

    try:
        threads = int(kv.get("THREADS", "8"))
    except ValueError:
        raise ConfError(f"THREADS must be an integer: {kv['THREADS']!r}")
    if threads < 1:
        raise ConfError(f"THREADS must be at least 1: {threads}")
    java_mem = kv.get("JAVA_MEM", "32G")
    if not (java_mem[:-1].isdigit() and java_mem[-1:].upper() in ("M", "G")):
        raise ConfError(f"JAVA_MEM format, e.g. 32G or 8000M → {java_mem!r}")

    res_dir = resolve_path(kv["GATK_RESOURCE_DIR"], base)
    resources = {k: resolve_path(kv[k], base) if kv.get(k) else res_dir / name
                 for k, name in GATK_RESOURCES.items()}

    beds = {"WGS": WGS_BED}
    for p, key in USER_BED_KEYS.items():
        beds[p] = resolve_path(kv[key], base) if kv.get(key) else None

    tools = {}
    for key, (env, default) in TOOL_KEYS.items():
        v = kv.get(key) or os.environ.get(env) or default
        tools[key] = str(resolve_path(v, base)) if ("/" in v or "\\" in v) else v   # resolve paths; bare names are looked up on PATH

    return Settings(
        path=conf_path,
        ref_fasta=resolve_path(kv["REF_FASTA"], base),
        gatk_resource_dir=res_dir,
        results_dir=resolve_path(kv["RESULTS_DIR"], base),
        pon_dir=resolve_path(kv["PON_DIR"], base) if kv.get("PON_DIR") else None,
        bam_dir=resolve_path(kv["BAM_DIR"], base) if kv.get("BAM_DIR") else None,
        threads=threads,
        java_mem=java_mem,
        beds=beds,
        resources=resources,
        tools=tools,
        ncbi_api_key=kv.get("NCBI_API_KEY") or os.environ.get("NCBI_API_KEY") or None,
        unknown_keys=sorted(set(kv) - KNOWN_KEYS),
    )
