# DeepPurityMap App

DeepPurityMap also comes as a browser app and gives the copy-number profile, gene-level calls and an HTML report for each sample

- **Install**: Docker and one launcher. The tools, the model and the databases come in one image.
- **Set up**: choose files in the app. Every value is checked, the sample sheet is filled from the BAM folder, and a missing PoN is built automatically.
- **Run**: one click runs every sample through all seven steps. Finished steps are skipped, and closing the browser does not stop the run.
- **Interpret**: read each sample's report in the app and ask an AI chat about it.

---

## 1. Prerequisite

### Software

**Docker** — the app, the analysis tools, the model and the databases all come in one Docker image.

| | Install first |
|---|---|
| **Linux (Ubuntu)** | [Docker Engine](https://docs.docker.com/engine/install/) |
| **Windows** | [Docker Desktop](https://docs.docker.com/desktop/setup/install/windows-install/) (WSL 2) |
| **macOS** | not yet tested |


A GPU is optional: without one, the purity model runs on the CPU.
To use an NVIDIA GPU on Linux, also install the [NVIDIA Container Toolkit](https://docs.nvidia.com/datacenter/cloud-native/container-toolkit/latest/install-guide.html); on Windows, Docker Desktop uses it as is.

### Data

| Item | Requirements | Where |
|---|---|---|
| **Reference FASTA (hg38)** | The same reference the BAMs were aligned to. Uncompressed `.fasta`, with `.fai` and `.dict` in the same folder | [Broad public bucket](https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/): `Homo_sapiens_assembly38.fasta` · `.fai` · `.dict` |
| **4 GATK resources** | In one folder, each with its `.tbi`:<br>`dbsnp_146.hg38.vcf.gz`<br>`Mills_and_1000G_gold_standard.indels.hg38.vcf.gz`<br>`Homo_sapiens_assembly38.known_indels.vcf.gz`<br>`1000G_omni2.5.hg38.vcf.gz` | dbSNP 146: [Broad FTP](ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/) or the [DDBJ mirror](https://ddbj.nig.ac.jp/public/public-human-genomes/GRCh38/fasta/)<br>Others: the Broad public bucket above |
| **BAM files** | hg38-aligned BAMs (the output of a standard alignment pipeline) | your data |
| **Capture kit BED** | WES and panel only: the kit maker's hg38 BED | the kit maker (WGS needs none) |
| **6+ normal BAMs** | Only for tumor-only samples: normals of the same capture kit, used to build the panel of normals (PoN) | your data |

The model and the annotation databases are already in the Docker image — no separate download is needed.

---

## 2. Get DeepPurityMap

Download **`DeepPurityMap-launcher.zip`** from the [Releases](https://github.com/jinoklee/DeepPurityMap/releases) page
and unpack it anywhere.

---

## 3. Running the analysis

### Start the app

Double-click the launcher in the unpacked `DeepPurityMap` folder:

| | Launcher |
|---|---|
| **Windows** | `DeepPurityMap.bat` |
| **Linux** | `DeepPurityMap.sh` (or run `./DeepPurityMap.sh` in a terminal) |

The app opens in your browser at `http://127.0.0.1:8800` (or the next free port).   
The first start downloads the DeepPurityMap image (about 7 GB; keep about 20 GB free on the Docker disk), later starts take a few seconds.  
The folders of your computer are available in the app, so use your usual paths (`D:\data\...` on Windows, `/data/...` on Linux).

### Closing the window and quitting the app

| | Close the window | Quit the app |
|---|---|---|
| How | close the browser tab or window | the **⏻** button at the top right of the app |
| The app | **keeps running** | stops |
| A running analysis | **keeps running** | stops. The app asks first. Run it again later and the finished steps are skipped |
| To come back | double-click the launcher again. It opens the running app | double-click the launcher. It starts the app again |

Quitting Docker also quits the app.

### Analysis tab

One project folder per analysis:

1. **Project folder** — create one, or open a recent project.
2. **Settings** — the reference, GATK resources, BED and BAM folder of section 1. Saving checks each value.
   Tumor-only samples need a PoN: give an existing PoN folder, or the folder of normal BAMs (`NORMAL_BAM_DIR`) —
   then the PoN is built automatically, first, when the analysis starts.
   With a BAM folder set, saving also fills the sample sheet (`sample_id`, BAM path, platform).
3. **Sample sheet** — complete `patient_id`, `tumor`/`normal` and cancer type, then save.
   A tumor with a normal of the same `patient_id` is analysed as paired, otherwise as tumor-only (needs a PoN).
4. **Check**, then **Run**. Progress is shown per sample and step.

The chat panel on the right takes you through the same steps. How it works depends on whether an LLM is set:

| | With an LLM | Without an LLM |
|---|---|---|
| Start | Set it in ⚙ at the top right (provider, model, API key), or type the model and key into the chat (e.g. `claude-sonnet-5-5 sk-ant-…`). The chat then takes over. | The chat starts a **step-by-step guide**: one question at a time (project folder → platform → reference → … → check → run). Each answer is checked and filled in on the left. |
| How you work | Just say it: *"Open D:\projects\lung"*, *"Download the reference to D:\ref\hg38"*, *"The BAMs are in D:\bam — 01A is the tumor, 11A the normal"*, *"Save it"*. | Answer in the input box or press the buttons; [Stop guide] lets you fill the left side yourself. |
| Check and run | Check and run are offered as cards and start only when you press [Apply]. After the check, the chat reports the result and offers the next step. | The guide runs the check and offers [Start analysis]. |

LLM choices: Anthropic Claude, OpenAI, Google Gemini, Qwen, or a local model (Ollama, LM Studio, vLLM). The API key is kept
only while the app runs.
**Data sent to a cloud LLM**: the settings, the sample sheet (sample IDs, BAM paths), folder listings and run status
(Analysis tab), and purity, copy-number segments and annotated genes (Results tab). BAM or other file contents are never
sent. Use pseudonymized sample IDs, or a local model (e.g. Ollama) to keep everything on your computer.

### Results tab

The samples of the open project are listed in a table. Click a sample to read its report.
With an LLM set, the chat panel answers questions about the report or the whole project (*"Summarize the purities"*,
*"Show chr8"*, *"Which genes are gained?"*); without one, the reports are read on their own.

---

## 4. Results

Each tumor sample gets its own folder in the project's `results/` folder.
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
| `4_purity/<sample_id>.purity.tsv` | **Predicted tumor purity** |
| `5_cnv/<sample_id>.segments_cn.tsv` | Purity-adjusted copy number per segment |
| `5_cnv/<sample_id>.genes_cn.tsv` | Copy number and gain/loss call per gene |
| `5_cnv/<sample_id>.cnv_summary.json` | Summary of the copy-number profile |
| `6_plots/<sample_id>_purity.png` | Purity figure (SNP depth–VAF distribution with the expected positions for each purity) |
| `6_plots/<sample_id>_genome.png` | Genome-wide copy-number figure |
| `6_plots/<sample_id>_chr<N>.png` | Copy-number figure per chromosome |
| `7_report/<sample_id>_report_en.html` | **The report** (`_report_ko.html`: Korean) |
| `7_report/<sample_id>_gene_classification.csv` | Altered genes with OncoKB · CIViC categories |
| `metadata.json` | Status, start/end time and duration of each step, purity, sample sheet values |
| `logs/<step>.log` · `logs/time.tsv` | One log per step and the time per step — look here when a step fails |

---

## 5. Viewing the results

**In the app** — Results tab (section 3).

**Or open the HTML report in a browser.** It is a single self-contained file. It contains a summary, the tumor purity with its SNP depth–VAF map, the copy-number profile per chromosome,
therapy database matches for altered genes, and related PubMed literature.

## Citation

To be added.

## License

DeepPurityMap is available **for academic and non-commercial research use only**.
Commercial use, including use in a commercial service or product, is not permitted without a separate license
(contact: <contact email>). DeepPurityMap is a research tool and must not be used for clinical decisions.
See [LICENSE](LICENSE) for the full terms.

The report uses third-party annotation data under their own terms:

- **OncoKB** gene list — free for academic research use only; commercial use requires a license from [OncoKB](https://www.oncokb.org/terms)
- **CIViC** therapies — [CC0](https://civicdb.org/about)
- **NCBI Gene** — public domain ([NCBI policies](https://www.ncbi.nlm.nih.gov/home/about/policies/))
