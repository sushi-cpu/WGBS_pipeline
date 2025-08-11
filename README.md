# Whole Genome Bisulphite Sequencing Pipelines

This repository contains bash and nextflow scripts for Whole Genome Bisulphite Sequencing (WGBS)
Analysis pipeline. 

---

## Features of wgbs.sh

### 1. Flexible Pipeline Modes
- **WGBS Mode**: Processes WGBS data with alignment, deduplication, and methylation extraction.
- **MeDIP-seq Mode**: Processes MeDIP-seq data with alignment and peak calling using MACS2 (`--medip`). [EXPERIMENTAL]
- **Specialized Modes**:
  - `--meth-report`: Generates a methylation summary CSV (CpG, CHG, CHH percentages, mean coverage) from existing WGBS results.
  - `--bedgraph-only`: Creates bedGraph and bigWig files from WGBS coverage files for visualization.
  - `--qc`: Runs FastQC on raw FASTQ files and exits.
  - `--tqc`: Performs trimming and post-trim FastQC, then exits.
  - `--resume`: Resumes processing from trimmed FASTQ files.

### 2. Input and Configuration
- **Flexible Input**:
  - Supports paired-end FASTQ files via `--r1` and `--r2` or auto-detected in the `data/` directory.
  - Validates FASTQ files for non-empty content and >1000 reads.
- **Customizable Parameters**:
  - `--adapter`: Adapter sequence (default: `AGATCGGAAGAGC`).
  - `--quality`: Quality threshold (default: 20).
  - `--min-length`: Minimum read length (default: 20).
  - `--threads`: Number of threads (default: 4), with SLURM auto-detection.
- **Reference Genome**:
  - Auto-locates FASTA files in the `genome/` directory.
  - Prepares Bismark (WGBS) or Bowtie2 (MeDIP-seq) indices if absent.

### 3. Quality Control and Preprocessing
- **Pre-Trim QC**: Runs FastQC on raw FASTQ files.
- **Trimming**: Uses Cutadapt to remove adapters and low-quality bases.
- **Post-Trim QC**: Runs FastQC on trimmed FASTQ files.

### 4. Alignment and Processing
- **WGBS Pipeline**:
  - Aligns reads using Bismark with Bowtie2.
  - Deduplicates PCR duplicates with `deduplicate_bismark`.
  - Extracts methylation calls (CpG, CHG, CHH) using `bismark_methylation_extractor`.
- **MeDIP-seq Pipeline**:
  - Aligns reads using Bowtie2.
  - Converts SAM to sorted/indexed BAM with Samtools.
  - Calls peaks with MACS2, producing `narrowPeak` and XLS files.

### 5. Output Generation
- **Directory Structure**: Organizes results in sample-specific `qc/`, `trimmed/`, `aligned/`, `methylation/` (WGBS), or `peaks/` (MeDIP-seq) directories.
- **Outputs**:
  - **WGBS**: Coverage files (`.cov`, `.cov.gz`), bedGraph, bigWig, and splitting reports.
  - **MeDIP-seq**: Sorted BAM, peak files, and MACS2 logs.
- **Reports**:
  - `methylation_report.csv`: Summarizes WGBS methylation metrics.
  - MultiQC HTML report for aggregated QC metrics.
- **Chromosome Sizes**: Generates a `.chrom.sizes` file for visualization.

### 6. Error Handling and Logging
- Uses `set -euo pipefail` for robust error handling.
- Logs detailed timestamps and messages to sample-specific log files and stdout.
- Validates input files and read counts.

### 7. Optimization and Flexibility
- `--force`: Overwrites existing outputs.
- `--clean`: Removes intermediate files to save space.
- SLURM integration for thread optimization.
- Multi-core support for FastQC, Cutadapt, Bismark, and Bowtie2.

### 8. Visualization Support
- Generates sorted bedGraph and bigWig files for WGBS data (if `bedGraphToBigWig` is available).

### 9. User-Friendly Interface
- `-h/--help`: Displays a detailed usage guide.
- Summarizes output locations and processing time per sample.

### 10. Portability and Extensibility
- Consistent directory structure (`data/`, `genome/`) for portability.
- Modular design for easy extension.
- Integrates standard tools (FastQC, Cutadapt, Bismark, Bowtie2, Samtools, MACS2, MultiQC).

---

## Requirements

- [Conda / Miniconda](https://docs.conda.io/en/latest/miniconda.html) installed on your system
- Bash shell (Linux/macOS)
- FASTQ files placed in the `data/` directory

--- 

## Installation

1. **Clone the repository**
   ```bash
   git clone https://github.com/sushi-cpu/WGBS_pipeline.git
   cd WGBS_pipeline

  2. **Create the Conda environment**
     ```bash
     conda env create -f environment.yml

  3. **Activate the environment**
     ```bash
     conda activate wgbs-core

---

## Before running the script

Make sure you have the following folder structure
```
WGBS_pipeline/
│
├── data/                     # Contains raw FASTQ files (paired-end or single-end)
│   ├── sample1_R1.fastq.gz
│   ├── sample1_R2.fastq.gz
│   ├── sample2_R1.fastq.gz
│   └── sample2_R2.fastq.gz
│
├── genome/                   # Contains the reference genome FASTA file
│   └── hg38.fa
│
├── script/                   # Contains the WGBS pipeline script
│   └── wgbs.sh
│   └── wgbs.nf
|
└── README.md                 # Documentation
```

## Usage
```bash
chmod +x wgbs.sh
```
```bash
./wgbs.sh
```

### Modularity  
```bash
./wgbs.sh --help
```
This shows 
```bash
Usage: ./wgbs.sh [options]
Options:
  --adapter                  Adapter sequence (default: AGATCGGAAGAGC)
  --quality                  Quality threshold (default: 20)
  --min-length               Minimum read length (default: 20)
  --threads                  Number of threads (default: 4)
  --force                    Overwrite existing output
  --clean                    Remove intermediate files
  --qc                       Run only FastQC (pre-trim), then exit
  --tqc                      Run trimming + post-trim FastQC, then exit
  --resume                   Resume pipeline from trimmed FASTQ files
  --r1                       Path to R1 FASTQ file
  --r2                       Path to R2 FASTQ file
  --meth-report              Generate methylation summary only (WGBS)
  --bedgraph-only            Generate bedGraph and bigWig files only (WGBS)
  --medip                    Run MeDIP-seq pipeline instead of WGBS
```

### Example usage
```bash
./wgbs.sh --r1 ../data/sample1_R1.fastq.gz --r2 ../data/sample1_R2.fastq.gz --quality 35
```
```bash
./wgbs.sh --medip #To trigger MeDIP-seq pipeline
```


## Features of wgbs.nf

- **Fully Modular Nextflow DSL2 Pipeline**  
  Designed for Whole Genome Bisulfite Sequencing (WGBS) analysis, implemented using **Nextflow DSL2** for reproducibility, scalability, and portability.

- **Automated Genome Preparation**  
  Automatically detects if the Bismark genome index already exists, skipping unnecessary re-generation.

- **Raw Read Quality Control (FastQC)**  
  Generates detailed QC reports for all raw FASTQ files before processing.

- **Read Trimming with Cutadapt**  
  Removes adapters, trims low-quality bases, and filters short reads based on configurable parameters:
  - `--adapter` sequence
  - `--quality` cutoff
  - `--minimum-length` filtering  
  Logs are saved for full reproducibility.

- **Post-Trim Quality Control**  
  Runs FastQC again on trimmed reads to ensure read quality improvement.

- **Bismark Alignment**  
  Aligns paired-end bisulfite reads against the prepared genome:
  - Multi-core support
  - Automatic BAM file renaming
  - Detailed alignment reports

- **Deduplication**  
  Removes PCR duplicates using `deduplicate_bismark` and outputs clean BAM files for downstream analysis.

- **Methylation Extraction**  
  Produces methylation reports in multiple formats:
  - Cytosine reports
  - CpG reports
  - BedGraph and coverage files
  - Splitting and M-bias reports

- **Customizable Parameters**  
  Modify adapters, quality thresholds, minimum length, threads, and genome path directly in the `nextflow.config` or via CLI.

- **Organized Output Structure**  
  Results are neatly arranged in per-sample folders:
```
results/
├── SAMPLE_ID/
    ├── 01_QC_reports
    ├── 02_trimmed
    ├── 03_trimmed_QC_reports
    ├── 04_aligned
    ├── 05_deduplicated
    └── 06_methylation
```
I have uploaded the sample output files in this repository under "nextflow results" folder 


### Software
Make sure the following tools are installed and accessible in your `$PATH` before running the pipeline:

- [Nextflow](https://www.nextflow.io/) (>= 22.10.0, DSL2 enabled)
- [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) (>= 0.11.9)
- [Cutadapt](https://cutadapt.readthedocs.io/) (>= 4.0)
- [Bismark](https://www.bioinformatics.babraham.ac.uk/projects/bismark/) (>= 0.24.0)  
  Requires either **Bowtie2** or **HISAT2** (Bowtie2 recommended: >= 2.4.5)
- [Samtools](http://www.htslib.org/) (>= 1.15)

### Hardware
- Minimum 4 CPU cores (more recommended for faster processing)
- Minimum 8 GB RAM (per sample for alignment)
- Sufficient disk space (WGBS data is large — plan for at least 2× the size of your raw FASTQ files)

### Genome Data
- A reference genome in FASTA format placed in the `genome/` directory
- The genome will be indexed by `bismark_genome_preparation` if not already prepared


### Usage 
It requires the same folder structure as wgbs.sh

```bash
nexflow run wgbs.nf
```

## 📜 Citation

If you use this pipeline in your research, please cite:

**WGBS Pipeline**  
> Sushant Nimbhorkar. *WGBS / MeDIP-seq Analysis Pipeline*. GitHub repository: [https://github.com/sushi-cpu/WGBS_pipeline](https://github.com/sushi-cpu/WGBS_pipeline), 2025.

**Dependencies**  
- Andrews S. *FastQC: A quality control tool for high throughput sequence data*.  
- Martin M. *Cutadapt removes adapter sequences from high-throughput sequencing reads*. EMBnet.journal, 17(1):10-12, 2011.  
- Krueger F, Andrews SR. *Bismark: a flexible aligner and methylation caller for Bisulfite-Seq applications*. Bioinformatics, 27(11):1571–1572, 2011.  
- Langmead B, Salzberg SL. *Fast gapped-read alignment with Bowtie 2*. Nat Methods. 9(4):357–359, 2012.  
- Li H, et al. *The Sequence Alignment/Map format and SAMtools*. Bioinformatics, 25(16):2078–2079, 2009.  
- Zhang Y, et al. *Model-based Analysis of ChIP-Seq (MACS)*. Genome Biology, 9(9):R137, 2008.  
- Ewels P, et al. *MultiQC: summarize analysis results for multiple tools and samples in a single report*. Bioinformatics, 32(19):3047–3048, 2016.  

---

## 📄 License

This project is licensed under the **MIT License** – you are free to use, modify, and distribute it with attribution.

---

## 👤 Author

**Sushant Nimbhorkar**  
📧 sushantnimbhorkar@gmail.com  
💻 [GitHub Profile](https://github.com/sushi-cpu)  
🏢 LinkedIn - [Sushant Nimbhorkar](https://www.linkedin.com/in/sushant-nimbhorkar-a57385242/)





