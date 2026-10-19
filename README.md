# DeepPurityMap

**Deep-learning tumor purity estimation from 2D SNP depth–VAF maps**

DeepPurityMap estimates tumor purity from sequencing data. SNP depth and variant allele fraction (VAF) across the genome are
drawn as a 2D map, and a deep-learning model reads that map to predict purity.

The same analysis also comes as a browser app, [Web App](https://github.com/jinoklee/DeepPurityMap/App/README.md):

---

## 1. Requirements

Linux with [conda](https://docs.conda.io/) — we recommend [Miniforge](https://github.com/conda-forge/miniforge).

### Get DeepPurityMap

```bash
git clone https://github.com/jinoklee/DeepPurityMap.git
cd DeepPurityMap
```

### Create the conda environments

Three environments are used: analysis tools and R, GATK, and the purity model.

```bash
# puritymap — samtools · bcftools · tabix · mosdepth · bedtools · R 4.3 · Python 3.10
conda env create -n puritymap -f environment.yml
conda install -n puritymap -c conda-forge openjdk=17 -y
conda run -n puritymap Rscript -e 'install.packages("PSCBS", repos="https://cloud.r-project.org")'

# puritymap-predict — PyTorch 2.6 · natten · transformers
conda create -n puritymap-predict python=3.10 -y
conda run -n puritymap-predict pip install -r requirements-predict.txt

# puritymap-gatk — GATK 4.6.2.0 and its Python environment
cd ~ && wget https://github.com/broadinstitute/gatk/releases/download/4.6.2.0/gatk-4.6.2.0.zip
unzip gatk-4.6.2.0.zip && cd gatk-4.6.2.0
conda env create -n puritymap-gatk -f gatkcondaenv.yml 

```

Check:

```bash
conda run -n puritymap Rscript -e 'library(PSCBS); library(vcfR); library(ggplot2); cat("R OK\n")'
conda run -n puritymap-predict python -c "import torch, natten, transformers; print('cuda:', torch.cuda.is_available())"
```

### Download the model

The model weights are attached to the [DeepPurityMap v1.0.0 release](https://github.com/jinoklee/DeepPurityMap/releases/tag/v1.0.0).
Download them into the DeepPurityMap folder (`resources/model/final_model_deploy.pth`):

```bash
cd ~/DeepPurityMap-CLI
wget -O resources/model/final_model_deploy.pth \
  https://github.com/jinoklee/DeepPurityMap/releases/download/v1.0.0/final_model_deploy.pth
```


---

## 2. Prerequisite

| Item | Requirements | Where to get it |
|---|---|---|
| **Reference FASTA (hg38)** | The same reference the BAMs were aligned to. Uncompressed `.fasta`, with `.fai` and `.dict` in the same folder | [Broad public bucket](https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/): `Homo_sapiens_assembly38.fasta` · `.fai` · `.dict` |
| **4 GATK resources** | In one folder, each with its `.tbi`:<br>`dbsnp_146.hg38.vcf.gz`<br>`Mills_and_1000G_gold_standard.indels.hg38.vcf.gz`<br>`Homo_sapiens_assembly38.known_indels.vcf.gz`<br>`1000G_omni2.5.hg38.vcf.gz` | dbSNP 146: [Broad FTP](ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/) or the [DDBJ mirror](https://ddbj.nig.ac.jp/public/public-human-genomes/GRCh38/fasta/)<br>Others: the Broad public bucket above |
| **BAM files** | hg38-aligned BAMs (the output of a standard alignment pipeline) | your data |
| **Capture kit BED** | WES and panel only: the kit maker's hg38 BED | the kit maker (WGS needs none) |
| **6+ normal BAMs** | Only for tumor-only samples: normals of the same capture kit, used to build the panel of normals (PoN) | your data |

Then make a **project folder** for the analysis and copy the two example files into it:

```bash
mkdir ~/project
cp ~/PurityMap/examples/puritymap.conf.example ~/project/puritymap.conf
cp ~/PurityMap/examples/samples.example.tsv    ~/project/samples.tsv
```

**`puritymap.conf`** — set the paths:

```ini
REF_FASTA=/data/ref/Homo_sapiens_assembly38.fasta
GATK_RESOURCE_DIR=/data/ref/gatk_bundle
RESULTS_DIR=results                       # relative to this file
WES_BED=/data/bed/capture_kit_hg38.bed    # PANEL_BED= for a panel; nothing for WGS
# PON_DIR=pon                             # tumor-only only
THREADS=8
JAVA_MEM=32G

# tool locations (see below)
GATK=/home/user/gatk-4.6.2.0/gatk
RSCRIPT=/home/user/miniforge3/envs/puritymap/bin/Rscript
PREDICT_PYTHON=/home/user/miniforge3/envs/puritymap-predict/bin/python
BEDTOOLS=/home/user/miniforge3/envs/puritymap/bin/bedtools
```

All four tool lines are required: the GATK environment has its own `python` and `Rscript`, which would be found
first otherwise. The locations differ per computer — print the lines for yours and paste them into `puritymap.conf`
(`GATK` is the `gatk` file in the folder where you unzipped GATK):

```bash
echo "GATK=$HOME/gatk-4.6.2.0/gatk"
echo "RSCRIPT=$(conda run -n puritymap which Rscript)"
echo "PREDICT_PYTHON=$(conda run -n puritymap-predict which python)"
echo "BEDTOOLS=$(conda run -n puritymap which bedtools)"
```

**`samples.tsv`** — one row per BAM, tab-separated:

```
sample_id	patient_id	type	bam	platform	cancer_type
P01_T	P01	tumor	/data/bam/P01_tumor.bam	WES	Stomach adenocarcinoma
P01_N	P01	normal	/data/bam/P01_normal.bam	WES
```

- `type`: `tumor` or `normal`. A tumor with a normal of the same `patient_id` is analysed as paired; otherwise as tumor-only.
- `platform`: `WES`, `Panel` or `WGS`.

---

## 3. Running the analysis

Activate the tool environments and use the `puritymap` environment's Python:

```bash
conda activate puritymap
conda activate --stack puritymap-gatk    # GATK's python tools need its environment first on PATH
cd ~/PurityMap
```

**1. Check** the settings, references, sheet, BAMs and tools. Fix what it reports and repeat until there are no errors.

```bash
python -m puritymap check --conf ~/project/puritymap.conf --sheet ~/project/samples.tsv
```

**2. Tumor-only samples only — build the PoN** once per capture kit, from a sheet that lists the normal BAMs:

```bash
pyrhon -m puritymap prepare-pon --conf ~/project/puritymap.conf --sheet ~/project/pon_normals.tsv --platform WES
```
**3. Run.** Run it in the background:

```bash
nohup pyrhon -m puritymap run --conf ~/project/puritymap.conf --sheet ~/project/samples.tsv > ~/project/run.log 2>&1 &
echo $! > ~/project/run.pid          # remember the run, to stop it later
tail -f ~/project/run.log            # follow the progress (Ctrl+C stops watching, not the run)
```

To stop the run, including the GATK and R processes it started:

```bash
kill -- -$(cat ~/project/run.pid)
```

- To resume a stopped or interrupted run, run the same command again — finished steps are skipped.
- Useful options: `--samples P01_T` (only these tumors), `--redo figure` (repeat a step), `--dry-run` (print the commands only).

---


## 4. Results

Each tumor sample gets its own folder under `RESULTS_DIR` (here `~/project/results/`):

The step folders are created in order as the run proceeds:

| Folder / file in `results/<sample_id>/` | Content |
|---|---|
| `1_genotype/final.vcf.gz` | Filtered germline SNPs (GATK HaplotypeCaller → variant filtering, PASS only) |
| `2_depth/<sample_id>.regions.bed.gz` | Read depth per region (mosdepth); for a paired sample also the matched normal's depth |
| `3_features/<sample_id>.feature.txt` | Normalized depth and VAF of each heterozygous SNP, with its segment |
| `3_features/<sample_id>.pscbs_seg.txt` | Copy-number segments |
| `3_features/<sample_id>.seg.png` | **2D SNP depth–VAF map** (segment level) — the image the model reads |
| `3_features/<sample_id>.color.png` | The same map at SNP level |
| `3_features/qc.tsv` | Quality-control values (normal reference used, heterozygous SNP count, noise) |
| `4_purity/<sample_id>.purity.tsv` | **Predicted tumor purity** (with the model name, version and device) |
| `4_purity/<sample_id>_purity.png` | **Purity figure** — SNP depth–VAF distribution with the expected positions for each purity |
| `metadata.json` | Status, start/end time and duration of each step, purity, sample sheet values |
| `logs/<step>.log` · `logs/time.tsv` | One log per step and the time per step — look here when a step fails |

---


## Citation

To be added.

## License

PurityMap is available **for academic and non-commercial research use only**.
Commercial use, including use in a commercial service or product, is not permitted without a separate license
(contact: <contact email>). PurityMap is a research tool and must not be used for clinical decisions.
See [LICENSE](LICENSE) for the full terms.
