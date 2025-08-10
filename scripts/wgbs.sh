#!/bin/bash
set -euo pipefail

# ────────────── Logging Helper ──────────────
log() {
  echo -e "[$(date +'%Y-%m-%d %H:%M:%S')] $*"
}

# ────────────── Configuration ──────────────
ADAPTER="AGATCGGAAGAGC"
QUALITY=20
MIN_LENGTH=20
THREADS=4
FORCE=0
CLEAN=0
QC_ONLY=0
TRIM_QC_ONLY=0
RESUME_FROM_TRIMMED=0
METH_REPORT_ONLY=0
BEDGRAPH_ONLY=0
MEDIP=0
R1_INPUT=""
R2_INPUT=""
BASE_DIR="$(cd "$(dirname "$0")/.." && pwd)"
RESULT_DIR_BASE="${BASE_DIR}"
DATA="${BASE_DIR}/data"
GENOME="${BASE_DIR}/genome"
GENOME_FASTA=$(find "$GENOME" -maxdepth 1 -name "*.fa" -o -name "*.fasta" | head -n 1)

if [[ -z "$GENOME_FASTA" ]]; then
  log "Could not find reference genome FASTA in $GENOME"
  exit 1
fi

# ────────────── Parse Arguments ──────────────
while [[ $# -gt 0 ]]; do
  case $1 in
    --adapter) ADAPTER="$2"; shift 2;;
    --quality) QUALITY="$2"; shift 2;;
    --min-length) MIN_LENGTH="$2"; shift 2;;
    --threads) THREADS="$2"; shift 2;;
    --force) FORCE=1; shift;;
    --clean) CLEAN=1; shift;;
    --qc) QC_ONLY=1; shift;;
    --tqc) TRIM_QC_ONLY=1; shift;;
    --resume) RESUME_FROM_TRIMMED=1; shift;;
    --r1) R1_INPUT="$2"; shift 2;;
    --r2) R2_INPUT="$2"; shift 2;;
    --meth-report) METH_REPORT_ONLY=1; shift;;
    --bedgraph-only) BEDGRAPH_ONLY=1; shift;;
    --medip) MEDIP=1; shift;;
    -h|--help)
      echo "Usage: $0 [options]"
      echo "Options:"
      echo "  --adapter STR              Adapter sequence (default: $ADAPTER)"
      echo "  --quality INT              Quality threshold (default: $QUALITY)"
      echo "  --min-length INT           Minimum read length (default: $MIN_LENGTH)"
      echo "  --threads INT              Number of threads (default: $THREADS)"
      echo "  --force                    Overwrite existing output"
      echo "  --clean                    Remove intermediate files"
      echo "  --qc                       Run only FastQC (pre-trim), then exit"
      echo "  --tqc                      Run trimming + post-trim FastQC, then exit"
      echo "  --resume                   Resume pipeline from trimmed FASTQ files"
      echo "  --r1 FILE                  Path to R1 FASTQ file"
      echo "  --r2 FILE                  Path to R2 FASTQ file"
      echo "  --meth-report              Generate methylation summary only (WGBS)"
      echo "  --bedgraph-only            Generate bedGraph and bigWig files only (WGBS)"
      echo "  --medip                    Run MeDIP-seq pipeline instead of WGBS"
      exit 0;;
    *) log "Unknown argument: $1"; exit 1;;
  esac
done

log "MEDIP flag set to: $MEDIP"

# ────────────── METH REPORT ONLY ──────────────
if [[ $METH_REPORT_ONLY -eq 1 ]]; then
  log "Generating methylation summary only"
  REPORT="${BASE_DIR}/methylation_report.csv"

  if [[ ! -f "$REPORT" ]]; then
    echo "Sample,CpG %,CHG %,CHH %,Mean CpG Coverage,M-bias QC,Report File" > "$REPORT"
  fi

  for SAMPLE_DIR in "${RESULT_DIR_BASE}"/*; do
    [[ -d "$SAMPLE_DIR" ]] || continue
    SAMPLE=$(basename "$SAMPLE_DIR")
    case "$SAMPLE" in
      scripts|genome|data|logs|*.sh|*.csv) continue ;;
    esac

    METH="${SAMPLE_DIR}/methylation"
    SPLIT_REPORT=$(find "$METH" -name '*splitting_report.txt' 2>/dev/null | head -n 1)
    COV_FILE=$(find "$METH" -type f \( -name '*.cov.gz' -o -name '*.cov' \) 2>/dev/null | head -n 1)
    MBIAS_FILE=$(find "$METH" -name '*.M-bias.txt' 2>/dev/null | head -n 1)

    if [[ ! -f "$SPLIT_REPORT" || ! -f "$COV_FILE" ]]; then
      log "Required files not found for $SAMPLE — skipping."
      continue
    fi

    CPG=$(grep "in CpG context" "$SPLIT_REPORT" | grep -o '[0-9.]\+%' | tr -d '%')
    CHG=$(grep "in CHG context" "$SPLIT_REPORT" | grep -o '[0-9.]\+%' | tr -d '%')
    CHH=$(grep "in CHH context" "$SPLIT_REPORT" | grep -o '[0-9.]\+%' | tr -d '%')

    if [[ "$COV_FILE" == *.gz ]]; then
      AVG_COV=$(gunzip -c "$COV_FILE" | awk '{ sum += $5 } END { if (NR > 0) printf("%.2f", sum / NR); else print "NA" }')
    else
      AVG_COV=$(awk '{ sum += $5 } END { if (NR > 0) printf("%.2f", sum / NR); else print "NA" }' "$COV_FILE")
    fi

    [[ -s "$MBIAS_FILE" ]] && MBIAS_FLAG="✅" || MBIAS_FLAG="Missing"

    TMP_REPORT="${REPORT}.tmp"
    grep -v "^$SAMPLE," "$REPORT" > "$TMP_REPORT"
    mv "$TMP_REPORT" "$REPORT"

    echo "$SAMPLE,$CPG,$CHG,$CHH,$AVG_COV,$MBIAS_FLAG,$SPLIT_REPORT" >> "$REPORT"
    log "Added methylation stats for $SAMPLE"
  done

  log "Methylation summary report written to $REPORT"
  exit 0
fi

# ────────────── BEDGRAPH ONLY ──────────────
if [[ $BEDGRAPH_ONLY -eq 1 ]]; then
  log "🧬 Generating bedGraph and bigWig files only"
  log "📂 Scanning sample directories in: $RESULT_DIR_BASE"

  for SAMPLE_DIR in "${RESULT_DIR_BASE}"/*; do
    [[ -d "$SAMPLE_DIR" ]] || continue
    SAMPLE=$(basename "$SAMPLE_DIR")
    case "$SAMPLE" in
      data|genome|multiqc_report|scripts|WGBS_docker) continue ;;
    esac

    METH_DIR="${SAMPLE_DIR}/methylation"
    [[ -d "$METH_DIR" ]] || continue

    COV_FILE=$(find "$METH_DIR" -type f \( -name "*.cov.gz" -o -name "*.cov" \) 2>/dev/null | head -n 1)
    [[ -f "$COV_FILE" ]] || continue

    TMP_COV=""
    if [[ "$COV_FILE" == *.gz ]]; then
      TMP_COV="${COV_FILE%.gz}"
      gunzip -c "$COV_FILE" > "$TMP_COV"
      COV_FILE="$TMP_COV"
    fi

    BEDGRAPH_FILE="${METH_DIR}/${SAMPLE}.CpG.bedGraph"
    BIGWIG_FILE="${METH_DIR}/${SAMPLE}.CpG.bw"

    awk 'BEGIN {OFS="\t"} {print $1, $2 - 1, $3, $4}' "$COV_FILE" > "$BEDGRAPH_FILE"
    log "📝 Generated bedGraph: $BEDGRAPH_FILE"

    sort -k1,1 -k2,2n "$BEDGRAPH_FILE" -o "$BEDGRAPH_FILE"
    log "🔃 Sorted bedGraph for $SAMPLE"

    if command -v bedGraphToBigWig &> /dev/null && [[ -f "$CHROM_SIZES" ]]; then
      bedGraphToBigWig "$BEDGRAPH_FILE" "$CHROM_SIZES" "$BIGWIG_FILE"
      log "📈 Created bigWig: $BIGWIG_FILE"
    fi

    [[ -n "$TMP_COV" && -f "$TMP_COV" ]] && rm -f "$TMP_COV"
  done

  log "Completed bedGraph-only processing"
  exit 0
fi

# ────────────── SLURM Threads Detection ──────────────
if [[ -n "${SLURM_JOB_ID:-}" ]]; then
  log "Running under SLURM (Job ID: $SLURM_JOB_ID)"
  THREADS=${SLURM_CPUS_PER_TASK:-$THREADS}
fi

# ────────────── Check and Prepare Genome ──────────────
if [[ $MEDIP -eq 0 ]]; then
  log "Checking Bismark genome preparation..."
  if [[ -f "$GENOME/Bisulfite_Genome/CT_conversion/genome_mfa.CT_conversion.fa" ]]; then
    log "Genome already prepared. Skipping genome preparation."
  else
    log "Preparing Bismark genome (this may take a while)..."
    bismark_genome_preparation --bowtie2 "$GENOME"
    log "Genome preparation complete."
  fi
else
  log "MeDIP-seq mode: Checking Bowtie2 index..."
  if ls "$GENOME"/genome*.bt2 >/dev/null 2>&1; then
    log "Bowtie2 index found. Skipping index preparation."
  else
    log "Preparing Bowtie2 index (this may take a while)..."
    bowtie2-build "$GENOME_FASTA" "$GENOME/genome"
    log "Bowtie2 index preparation complete."
  fi
fi

# ────────────── Find Sample Pairs ──────────────
if [[ -n "$R1_INPUT" && -n "$R2_INPUT" ]]; then
  R1_FILES=($R1_INPUT)
  log "Processing specified FASTQ pair: $R1_INPUT, $R2_INPUT"
else
  log "Scanning for FASTQ pairs in $DATA"
  R1_FILES=($(find "$DATA" -type f \( -name '*_R1*.fastq*' -o -name '*_1*.fastq*' \) | sort))
fi

if [[ ${#R1_FILES[@]} -eq 0 ]]; then
  log "No FASTQ R1 files found."
  exit 1
fi

# ────────────── Process Each Sample ──────────────
for R1 in "${R1_FILES[@]}"; do
  if [[ -n "$R2_INPUT" ]]; then
    R2="$R2_INPUT"
  else
    R2="${R1/_R1/_R2}"
    R2="${R2/_1/_2}"
  fi

  if [[ ! -f "$R2" ]]; then
    log "Skipping: Pair not found for $R1"
    continue
  fi

  # Validate FASTQ files
  if [[ ! -s "$R1" || ! -s "$R2" ]]; then
    log "Skipping: Empty or invalid FASTQ files: $R1, $R2"
    continue
  fi

  # Check read count for FASTQ files
  if [[ "$R1" == *.gz ]]; then
    if ! gunzip -c "$R1" >/dev/null 2>&1; then
      log "Error: Cannot decompress $R1. Skipping."
      continue
    fi
    READ_COUNT=$(gunzip -c "$R1" | wc -l)
  else
    READ_COUNT=$(wc -l < "$R1")
  fi
  READ_COUNT=$((READ_COUNT / 4))
  if [[ $READ_COUNT -lt 1000 ]]; then
    log "Warning: Low read count ($READ_COUNT) in $R1. Expected >1000 reads for MeDIP-seq/WGBS."
  fi

  SAMPLE=$(basename "$R1" | sed -E 's/(_R1|_1)\.fastq(.gz)?$//')
  log "Checking sample: $SAMPLE"

  SAMPLE_DIR="${BASE_DIR}/${SAMPLE}"
  QC="${SAMPLE_DIR}/qc"
  TRIMMED="${SAMPLE_DIR}/trimmed"
  ALIGNED="${SAMPLE_DIR}/aligned"
  LOG_FILE="${SAMPLE_DIR}/pipeline.log"

  # Check if sample is already processed
  if [[ $MEDIP -eq 1 ]]; then
    PEAKS="${SAMPLE_DIR}/peaks"
    mkdir -p "$QC/post_trim" "$TRIMMED" "$ALIGNED" "$PEAKS"
    NARROW_PEAK="${PEAKS}/${SAMPLE}_peaks.narrowPeak"
    PEAKS_XLS="${PEAKS}/${SAMPLE}_peaks.xls"
    MACS_LOG="${PEAKS}/${SAMPLE}_macs2.log"
    if [[ -f "$NARROW_PEAK" && -f "$PEAKS_XLS" && -f "$MACS_LOG" && $FORCE -eq 0 ]]; then
      log "Sample $SAMPLE already processed (peaks, XLS, and log found). Skipping."
      continue
    elif [[ -f "$NARROW_PEAK" || -f "$PEAKS_XLS" || -f "$MACS_LOG" ]]; then
      log "Partial MeDIP-seq outputs found for $SAMPLE:"
      [[ -f "$NARROW_PEAK" ]] && log "  • $NARROW_PEAK"
      [[ -f "$PEAKS_XLS" ]] && log "  • $PEAKS_XLS"
      [[ -f "$MACS_LOG" ]] && log "  • $MACS_LOG"
      if [[ $FORCE -eq 0 ]]; then
        log "Use --force to overwrite partial outputs."
        continue
      fi
    fi
  else
    METH="${SAMPLE_DIR}/methylation"
    mkdir -p "$QC/post_trim" "$TRIMMED" "$ALIGNED" "$METH"
    COV_FILE=$(find "$METH" -type f \( -name "*.cov.gz" -o -name "*.cov" \) 2>/dev/null | head -n 1)
    if [[ -f "$COV_FILE" && $FORCE -eq 0 ]]; then
      log "Sample $SAMPLE already processed (methylation data found). Skipping."
      continue
    fi
  fi

  log "🚀 Starting pipeline for sample: $SAMPLE"
  START_TIME=$(date +%s)

  if [[ $MEDIP -eq 0 && $RESUME_FROM_TRIMMED -eq 1 ]]; then
    if [[ ! -f "${TRIMMED}/${SAMPLE}_R1.trimmed.fastq.gz" || ! -f "${TRIMMED}/${SAMPLE}_R2.trimmed.fastq.gz" ]]; then
      log "Trimmed files not found for $SAMPLE. Use --force to generate new trimmed files."
      continue
    fi
    log "⏩ Resuming from trimmed files. Skipping pre-trim QC and trimming."
  else
    log "1. FastQC (Pre-trim)"
    fastqc --threads "$THREADS" "$R1" "$R2" -o "$QC" 2>> "$LOG_FILE"

    if [[ $QC_ONLY -eq 1 ]]; then
      log "QC-only mode enabled. Skipping rest of pipeline for $SAMPLE."
      echo "QC report available at: $QC"
      continue
    fi

    if [[ $MEDIP -eq 0 && -f "${TRIMMED}/${SAMPLE}_R1.trimmed.fastq.gz" && -f "${TRIMMED}/${SAMPLE}_R2.trimmed.fastq.gz" && $FORCE -eq 0 ]]; then
      log "Trimmed files already exist for $SAMPLE. Use --force to overwrite or --resume-from-trimmed to proceed."
      continue
    fi

    log "2. Trimming with Cutadapt"
    cutadapt \
        -j "$THREADS" \
        -a "$ADAPTER" -A "$ADAPTER" \
        -q "$QUALITY,$QUALITY" \
        --minimum-length "$MIN_LENGTH" \
        -o "${TRIMMED}/${SAMPLE}_R1.trimmed.fastq.gz" \
        -p "${TRIMMED}/${SAMPLE}_R2.trimmed.fastq.gz" \
        "$R1" "$R2" > "${TRIMMED}/${SAMPLE}_cutadapt.log" 2>&1

    log "3. FastQC (Post-trim)"
    fastqc --threads "$THREADS" "${TRIMMED}/${SAMPLE}_R1.trimmed.fastq.gz" "${TRIMMED}/${SAMPLE}_R2.trimmed.fastq.gz" -o "$QC/post_trim" 2>> "$LOG_FILE"

    if [[ $TRIM_QC_ONLY -eq 1 ]]; then
      log "Trim+QC-only mode enabled. Skipping rest of pipeline for $SAMPLE."
      echo "Post-trim QC report available at: $QC/post_trim"
      continue
    fi
  fi

  if [[ $MEDIP -eq 1 ]]; then
    log "4. Alignment with Bowtie2 (MeDIP-seq)"
    bowtie2 \
      -p "$THREADS" \
      -x "$GENOME/genome" \
      -1 "${TRIMMED}/${SAMPLE}_R1.trimmed.fastq.gz" \
      -2 "${TRIMMED}/${SAMPLE}_R2.trimmed.fastq.gz" \
      -S "${ALIGNED}/${SAMPLE}_aligned.sam" > "${ALIGNED}/${SAMPLE}_bowtie2.log" 2>&1

    log "5. Convert SAM to BAM and sort"
    samtools view -bS "${ALIGNED}/${SAMPLE}_aligned.sam" > "${ALIGNED}/${SAMPLE}_aligned.bam" 2>> "$LOG_FILE"
    samtools sort "${ALIGNED}/${SAMPLE}_aligned.bam" -o "${ALIGNED}/${SAMPLE}_aligned_sorted.bam" 2>> "$LOG_FILE"
    samtools index "${ALIGNED}/${SAMPLE}_aligned_sorted.bam" 2>> "$LOG_FILE"

    log "6. Peak Calling with MACS2"
    macs2 callpeak \
      -t "${ALIGNED}/${SAMPLE}_aligned_sorted.bam" \
      -f BAMPE \
      -g hs \
      --outdir "$PEAKS" \
      -n "$SAMPLE" \
      --nomodel \
      --extsize 200 \
      -q 0.01 > "${PEAKS}/${SAMPLE}_macs2.log" 2>&1

    if [[ $CLEAN -eq 1 ]]; then
      log "🧹 Cleaning up intermediate files"
      rm -f "${TRIMMED}/${SAMPLE}_R1.trimmed.fastq.gz" "${TRIMMED}/${SAMPLE}_R2.trimmed.fastq.gz"
      rm -f "${ALIGNED}/${SAMPLE}_aligned.sam" "${ALIGNED}/${SAMPLE}_aligned.bam"
    fi

    END_TIME=$(date +%s)
    TOTAL_TIME=$((END_TIME - START_TIME))

    log "MeDIP-seq pipeline complete!"
    echo "Outputs generated for sample $SAMPLE:"
    echo "• Sample directory: $SAMPLE_DIR"
    echo "• Quality control: $QC"
    echo "• Trimmed reads: $TRIMMED"
    echo "• Alignment files: $ALIGNED"
    echo "• Peak files: $PEAKS"
    echo "• Log file: $LOG_FILE"
    echo -n "• Total time taken: "
    printf '%02dh:%02dm:%02ds\n' $((TOTAL_TIME/3600)) $(((TOTAL_TIME/60)%60)) $((TOTAL_TIME%60))

    continue
  fi

  log "4. Alignment with Bismark"
  bismark \
    --multicore "$((THREADS/2))" \
    --genome "$GENOME" \
    -1 "${TRIMMED}/${SAMPLE}_R1.trimmed.fastq.gz" \
    -2 "${TRIMMED}/${SAMPLE}_R2.trimmed.fastq.gz" \
    -o "$ALIGNED" 2>> "$LOG_FILE"

  BAM="${ALIGNED}/${SAMPLE}_R1.trimmed_bismark_bt2_pe.bam"

  log "5. Deduplication"
  cd "$ALIGNED"
  deduplicate_bismark -p --bam "$(basename "$BAM")" 2>> "$LOG_FILE"
  DEDUP="${BAM%.bam}.deduplicated.bam"
  cd "$BASE_DIR"

  log "6. Methylation Extraction"
  bismark_methylation_extractor \
    --paired-end --gzip --bedGraph --cytosine_report \
    --genome_folder "$GENOME" \
    "$DEDUP" \
    -o "$METH" 2>> "$LOG_FILE"

  if [[ $CLEAN -eq 1 ]]; then
    log "🧹 Cleaning up trimmed files"
    rm -f "${TRIMMED}/${SAMPLE}_R1.trimmed.fastq.gz" "${TRIMMED}/${SAMPLE}_R2.trimmed.fastq.gz"
  fi

  END_TIME=$(date +%s)
  TOTAL_TIME=$((END_TIME - START_TIME))

  log "WGBS pipeline complete!"
  echo "Outputs generated for sample $SAMPLE:"
  echo "• Sample directory: $SAMPLE_DIR"
  echo "• Quality control: $QC"
  echo "• Trimmed reads: $TRIMMED"
  echo "• Alignment files: $ALIGNED"
  echo "• Methylation data: $METH"
  echo "• Log file: $LOG_FILE"
  echo -n "• Total time taken: "
  printf '%02dh:%02dm:%02ds\n' $((TOTAL_TIME/3600)) $(((TOTAL_TIME/60)%60)) $((TOTAL_TIME%60))
done

# ────────────── MultiQC Summary Report ──────────────
log "Generating MultiQC report"
MULTIQC_DIR="${BASE_DIR}/multiqc_report"
mkdir -p "$MULTIQC_DIR"
multiqc "$BASE_DIR" -o "$MULTIQC_DIR" --force 2>> "$BASE_DIR/multiqc.log"
log "MultiQC report generated at: $MULTIQC_DIR/multiqc_report.html"

# ────────────── WGBS Methylation Report ──────────────
if [[ $MEDIP -eq 0 ]]; then
  log "Generating Methylation Report"
  REPORT="${BASE_DIR}/methylation_report.csv"

  if [[ ! -f "$REPORT" ]]; then
    echo "Sample,CpG %,CHG %,CHH %,Mean CpG Coverage,M-bias QC,Report File" > "$REPORT"
  fi

  for SAMPLE_DIR in "${RESULT_DIR_BASE}"/*; do
    [[ -d "$SAMPLE_DIR" ]] || continue
    SAMPLE=$(basename "$SAMPLE_DIR")
    case "$SAMPLE" in
      scripts|genome|data|logs|*.sh|*.csv) continue ;;
    esac

    METH="${SAMPLE_DIR}/methylation"
    [[ -d "$METH" ]] || continue

    SPLIT_REPORT=$(find "$METH" -name '*splitting_report.txt' 2>/dev/null | head -n 1)
    COV_FILE=$(find "$METH" -type f \( -name '*.cov.gz' -o -name '*.cov' \) 2>/dev/null | head -n 1)
    MBIAS_FILE=$(find "$METH" -name '*.M-bias.txt' 2>/dev/null | head -n 1)

    if [[ ! -f "$SPLIT_REPORT" || ! -f "$COV_FILE" ]]; then
      log "Required files not found for sample $SAMPLE. Skipping report entry."
      continue
    fi

    CPG=$(grep "in CpG context" "$SPLIT_REPORT" | grep -o '[0-9.]\+%' | tr -d '%')
    CHG=$(grep "in CHG context" "$SPLIT_REPORT" | grep -o '[0-9.]\+%' | tr -d '%')
    CHH=$(grep "in CHH context" "$SPLIT_REPORT" | grep -o '[0-9.]\+%' | tr -d '%')

    if [[ "$COV_FILE" == *.gz ]]; then
      if gunzip -c "$COV_FILE" &>/dev/null; then
        AVG_COV=$(gunzip -c "$COV_FILE" | awk '{ sum += $5 } END { if (NR > 0) printf("%.2f", sum / NR); else print "NA" }')
      else
        log "Could not decompress $COV_FILE. Skipping entry."
        continue
      fi
    else
      AVG_COV=$(awk '{ sum += $5 } END { if (NR > 0) printf("%.2f", sum / NR); else print "NA" }' "$COV_FILE")
    fi

    [[ -s "$MBIAS_FILE" ]] && MBIAS_FLAG="" || MBIAS_FLAG=" Missing"

    TMP_REPORT="${REPORT}.tmp"
    grep -v "^$SAMPLE," "$REPORT" > "$TMP_REPORT"
    mv "$TMP_REPORT" "$REPORT"

    echo "$SAMPLE,$CPG,$CHG,$CHH,$AVG_COV,$MBIAS_FLAG,$SPLIT_REPORT" >> "$REPORT"
  done

  log "Methylation report generated at: $REPORT"
fi

# ────────────── Genome Chrom Sizes ──────────────
CHROM_SIZES="${GENOME_FASTA}.chrom.sizes"
if [[ ! -f "$CHROM_SIZES" ]]; then
  log "Generating chromosome sizes from $GENOME_FASTA"
  samtools faidx "$GENOME_FASTA"
  cut -f1,2 "${GENOME_FASTA}.fai" > "$CHROM_SIZES"
  log "Chrom sizes written to $CHROM_SIZES"
else
  log "Chrom sizes file found: $CHROM_SIZES"
fi

# ────────────── Generate bedGraph and bigWig (WGBS only) ──────────────
if [[ $MEDIP -eq 0 ]]; then
  log "Generating bedGraph and bigWig tracks from .cov files"
  find "$BASE_DIR" -mindepth 1 -maxdepth 1 -type d | while read -r SAMPLE_DIR; do
    METH_DIR="${SAMPLE_DIR}/methylation"
    SAMPLE=$(basename "$SAMPLE_DIR")
    if [[ ! -d "$METH_DIR" ]]; then
      continue
    fi

    COV_FILE=$(find "$METH_DIR" -type f \( -name "*.cov.gz" -o -name "*.cov" \) 2>/dev/null | head -n 1)
    if [[ ! -f "$COV_FILE" ]]; then
      log "No .cov file found for $SAMPLE. Skipping."
      continue
    fi

    TMP_COV=""
    if [[ "$COV_FILE" == *.gz ]]; then
      TMP_COV="${COV_FILE%.gz}"
      gunzip -c "$COV_FILE" > "$TMP_COV"
      COV_FILE="$TMP_COV"
    fi

    BEDGRAPH_FILE="${METH_DIR}/${SAMPLE}.CpG.bedGraph"
    BIGWIG_FILE="${METH_DIR}/${SAMPLE}.CpG.bw"

    awk 'BEGIN{OFS="\t"} {print $1, $2 - 1, $3, $4}' "$COV_FILE" > "$BEDGRAPH_FILE"
    log "Created bedGraph: $BEDGRAPH_FILE"

    sort -k1,1 -k2,2n "$BEDGRAPH_FILE" -o "$BEDGRAPH_FILE"

    if command -v bedGraphToBigWig &> /dev/null; then
      bedGraphToBigWig "$BEDGRAPH_FILE" "$CHROM_SIZES" "$BIGWIG_FILE"
      log "Created bigWig: $BIGWIG_FILE"
    else
      log "Skipping bigWig generation — bedGraphToBigWig not found"
    fi

    [[ -n "$TMP_COV" && -f "$TMP_COV" ]] && rm -f "$TMP_COV"
  done
fi

log " All samples processed."