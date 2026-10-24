#!/usr/bin/env bash
# =============================================================================
# 1_genotype_depth.sh   Track 1 steps 1 genotype · 2 depth (one sample)
#
#   tumor BAM → GATK HC → SNP → NVScoreVariants → FilterVariantTranches → PASS
#             → depth regions → mosdepth (tumor [+ matched normal])
#
#   Usage
#     1_genotype_depth.sh -t tumor.bam -s SAMPLE_ID -b intervals.bed -o RESULTS_DIR/SAMPLE_ID \
#                       [-n normal.bam -N NORMAL_ID] [-d depth_regions.bed]
#
#   Depth regions: by default the tumor's PASS SNP positions ±50bp (WES/Panel).
#   With -d (WGS) a fixed BED is used, so tumor and normal depth are on exactly the same windows.
#   Environment variables (filled from puritymap.conf by python -m puritymap)
#     REF_FASTA                         required
#     DBSNP MILLS KNOWN_INDELS OMNI     default: Broad standard file names inside GATK_RESOURCE_DIR
#     THREADS (8)  JAVA_MEM (32G)  GATK (gatk)
#
#   Output (under -o)
#     1_genotype/  genotype · snp · filtered · filtered_n1d · final.vcf.gz (+ indexes)
#     2_depth/     [snv.100bp.bed], {SAMPLE_ID}.regions.bed.gz, [{NORMAL_ID}.regions.bed.gz]
#     logs/        time.tsv
#
#   Finished steps are skipped (final.vcf.gz.tbi · *.regions.bed.gz.csi mark completion).
#   To rerun a step, delete its folder.
# =============================================================================
set -Eeuo pipefail                         # -E: failures inside functions also reach the ERR trap

usage() {
    sed -n '8,13p' "$0" | sed 's/^#//'
    exit "${1:-1}"
}
die() { echo "ERROR: $*" >&2; exit 1; }
now() { date '+%Y-%m-%d %H:%M:%S'; }
log() { echo "[$(date '+%H:%M:%S')] $*"; }

# ── Arguments ────────────────────────────────────────────────────────────────
tumor="" sample="" bed="" out="" normal="" normal_id="" depth_bed=""
while getopts "t:s:b:o:n:N:d:h" opt; do
    case ${opt} in
        t) tumor=${OPTARG} ;;
        s) sample=${OPTARG} ;;
        b) bed=${OPTARG} ;;
        o) out=${OPTARG} ;;
        n) normal=${OPTARG} ;;
        N) normal_id=${OPTARG} ;;
        d) depth_bed=${OPTARG} ;;
        h) usage 0 ;;
        *) usage ;;
    esac
done

[[ -n "${tumor}" && -n "${sample}" && -n "${bed}" && -n "${out}" ]] || { echo "ERROR: -t -s -b -o are required" >&2; usage; }
[[ -f "${tumor}" ]] || die "tumor BAM not found: ${tumor}"
[[ -f "${bed}" ]]   || die "intervals BED not found: ${bed}"
if [[ -n "${normal}" ]]; then
    [[ -f "${normal}" ]]   || die "normal BAM not found: ${normal}"
    [[ -n "${normal_id}" ]] || die "-n requires -N NORMAL_ID"
fi
[[ -z "${depth_bed}" || -f "${depth_bed}" ]] || die "depth regions BED not found: ${depth_bed}"

# ── Settings ─────────────────────────────────────────────────────────────────
: "${REF_FASTA:?REF_FASTA environment variable is required}"
[[ -f "${REF_FASTA}" ]] || die "REF_FASTA not found: ${REF_FASTA}"
res="${GATK_RESOURCE_DIR:-}"
DBSNP=${DBSNP:-${res}/dbsnp_146.hg38.vcf.gz}
MILLS=${MILLS:-${res}/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz}
KNOWN_INDELS=${KNOWN_INDELS:-${res}/Homo_sapiens_assembly38.known_indels.vcf.gz}
OMNI=${OMNI:-${res}/1000G_omni2.5.hg38.vcf.gz}
for f in "${DBSNP}" "${MILLS}" "${KNOWN_INDELS}" "${OMNI}"; do
    [[ -f "${f}" ]] || die "GATK resource not found: ${f}  (check GATK_RESOURCE_DIR or the individual variables)"
done
THREADS=${THREADS:-8}
JAVA_MEM=${JAVA_MEM:-32G}
GATK=${GATK:-gatk}
for t in "${GATK}" samtools bcftools tabix mosdepth; do
    command -v "${t}" >/dev/null 2>&1 || die "tool not found: ${t}"
done

gdir="${out}/1_genotype"
ddir="${out}/2_depth"
mkdir -p "${gdir}" "${ddir}" "${out}/logs"
TIME="${out}/logs/time.tsv"
[[ -f "${TIME}" ]] || echo -e "sample_id\tstep\tstart\tend\telapsed_s" > "${TIME}"
timed() {                                  # timed <step name> <command...>
    local step=$1; shift
    local s0 t0; s0=$(now); t0=$(date +%s)
    log "${step} start"
    "$@"
    echo -e "${sample}\t${step}\t${s0}\t$(now)\t$(( $(date +%s) - t0 ))" >> "${TIME}"
    log "${step} done ($(( $(date +%s) - t0 ))s)"
}
trap 'echo "ERROR: ${sample} — ${BASH_COMMAND} (line ${LINENO}) failed" >&2' ERR

log "START ${sample}  tumor=$(basename "${tumor}")${normal:+  normal=$(basename "${normal}")}"

# ── BAM index ────────────────────────────────────────────────────────────────
index_bam() {
    local bam=$1
    if [[ ! -f "${bam}.bai" && ! -f "${bam%.bam}.bai" ]]; then
        log "creating index: $(basename "${bam}")"
        samtools index -@ "${THREADS}" "${bam}"
    fi
}
index_bam "${tumor}"
[[ -n "${normal}" ]] && index_bam "${normal}"

# ── 1 genotype ───────────────────────────────────────────────────────────────
genotype() {
    local g=${gdir}
    "${GATK}" --java-options "-Xmx${JAVA_MEM}" HaplotypeCaller \
        --reference "${REF_FASTA}" \
        --input "${tumor}" \
        --output "${g}/genotype.vcf.gz" \
        --intervals "${bed}" \
        --native-pair-hmm-threads "${THREADS}" \
        -G StandardAnnotation \
        -G StandardHCAnnotation \
        -GQB 10 -GQB 20 -GQB 30 -GQB 60

    "${GATK}" SelectVariants -R "${REF_FASTA}" -V "${g}/genotype.vcf.gz" \
        --select-type-to-include SNP \
        -O "${g}/snp.vcf.gz"

    "${GATK}" NVScoreVariants \
        -R "${REF_FASTA}" \
        -V "${g}/snp.vcf.gz" \
        -O "${g}/filtered.vcf.gz"
    "${GATK}" IndexFeatureFile -I "${g}/filtered.vcf.gz"

    "${GATK}" FilterVariantTranches \
        -V "${g}/filtered.vcf.gz" \
        --resource "${DBSNP}" \
        --resource "${OMNI}" \
        --output "${g}/filtered_n1d.vcf.gz" \
        --info-key CNN_1D

    # PASS only. final.vcf.gz.tbi appears only if everything succeeded (= completion marker)
    bcftools view -f PASS "${g}/filtered_n1d.vcf.gz" -Oz -o "${g}/final.vcf.gz.tmp"
    mv "${g}/final.vcf.gz.tmp" "${g}/final.vcf.gz"
    tabix -f -p vcf "${g}/final.vcf.gz"
}
if [[ -f "${gdir}/final.vcf.gz.tbi" ]]; then
    log "genotype skipped (exists: 1_genotype/final.vcf.gz)"
else
    timed genotype genotype
fi

# ── 2 depth ──────────────────────────────────────────────────────────────────
if [[ -n "${depth_bed}" ]]; then
    # Fixed regions (WGS): tumor and normal on the same windows
    by_bed="${depth_bed}"
else
    # SNP position ±50bp regions (start clipped at 0 so it never goes negative)
    by_bed="${ddir}/snv.100bp.bed"
    if [[ ! -s "${by_bed}" || "${by_bed}" -ot "${gdir}/final.vcf.gz" ]]; then
        bcftools view -v snps "${gdir}/final.vcf.gz" \
            | awk 'BEGIN{OFS="\t"} !/^#/ { s = $2 - 50; if (s < 0) s = 0; print $1, s, $2 + 50 }' \
            > "${by_bed}.tmp"
        mv "${by_bed}.tmp" "${by_bed}"
    fi
    [[ -s "${by_bed}" ]] || die "0 PASS SNPs → check 1_genotype/final.vcf.gz"
fi
log "depth regions: ${by_bed}"

depth_one() {                              # depth_one <id> <bam>
    local id=$1 bam=$2
    if [[ -f "${ddir}/${id}.regions.bed.gz.csi" ]]; then
        log "depth skipped (exists: 2_depth/${id}.regions.bed.gz)"
        return
    fi
    timed "depth:${id}" mosdepth --by "${by_bed}" --threads "${THREADS}" -n "${ddir}/${id}" "${bam}"
}
depth_one "${sample}" "${tumor}"
[[ -n "${normal}" ]] && depth_one "${normal_id}" "${normal}"

log "DONE ${sample}"
