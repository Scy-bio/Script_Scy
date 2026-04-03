#!/bin/bash
# A dedicated Bash script for GATK RNA-seq SNP calling pipeline.
# It assumes STAR alignment has already been performed.

# Exit immediately if any command exits with a non-zero status.
set -e

echo "--- Starting GATK RNA-seq Variant Calling Pipeline ---"
echo "Current working directory: $(pwd)"
echo ""

# --- Default Configuration Variables ---
# These can be overridden by command-line arguments or adjusted manually.

# Input file paths (will be set by command-line arguments)
GENOME_FASTA=""   # Path to the genome FASTA file (e.g., region.fa)
SAMPLE_LIST=""    # Path to the sample list file (contains paths to STAR BAMs)

# General CPU Threads
SAMTOOLS_VIEW_THREADS=20  # Threads for samtools view operations

# GATK Java Options (adjust -Xmx for memory and -XX:ParallelGCThreads for cores)
# Ensure -Djava.io.tmpdir points to a directory with sufficient space for temporary files.
# These values are examples; adjust them based on your system's resources.
MARKDUP_JAVA_OPTIONS="-Xmx12g -XX:ParallelGCThreads=20 -Djava.io.tmpdir=${TEMP_GATK_DIR}" # Example: 8GB RAM, 20 threads for MarkDuplicates
SPLITNCIGAR_JAVA_OPTIONS="-Xmx12g -XX:ParallelGCThreads=10 -Djava.io.tmpdir=${TEMP_GATK_DIR}" # Example: 4GB RAM, 10 threads for SplitNCigarReads
HAPLOTYPECALLER_JAVA_OPTIONS="-Xmx12g -Djava.io.tmpdir=${TEMP_GATK_DIR}" # Example: 8GB RAM for HaplotypeCaller

# --- Output Directory Structure ---
# All generated output will be stored within this base directory.
GATK_OUTPUT_BASE="test_result/01.RNA-seq_variant_result"
REF_GENOME_DIR="${GATK_OUTPUT_BASE}/00.ref_genome" # New: Stores prepared reference genome files
BAM_PROCESSING_DIR="${GATK_OUTPUT_BASE}/01.processed_bams" # Stores uniquely mapped, RG-added, marked-dup BAMs
SPLIT_BAM_DIR="${GATK_OUTPUT_BASE}/02.split_bams"        # Stores SplitNCigarReads output BAMs
GVCF_DIR="${GATK_OUTPUT_BASE}/03.gvcf_calls"              # Stores GVCF files
TEMP_GATK_DIR="${GATK_OUTPUT_BASE}/tmp_gatk"              # Temporary directory for GATK operations

# --- Usage Function ---
usage() {
  echo "Usage: $0 --genome <path/to/region.fa> --rnaseq <path/to/list_file> [--help]"
  echo ""
  echo "Options:"
  echo "  --genome    Path to the genome FASTA file (e.g., data/region.fa)."
  echo "  --rnaseq    Path to the sample list file. Each line should be:"
  echo "              SampleID<tab>Path/To/SampleID.Aligned.sortedByCoord.out.bam"
  echo "  --help      Display this help message and exit."
  echo ""
  echo "Example:"
  echo "  ./run_gatk_pipeline.sh --genome data/region.fa --rnaseq list_of_aligned_bams.txt"
  exit 1
}

# --- Parse Command Line Arguments ---
OPTIONS=$(getopt -o "" -l genome:,rnaseq:,help -- "$@")

if [ $? -ne 0 ]; then
  echo "Error: Failed to parse command-line arguments." >&2
  usage
fi

eval set -- "$OPTIONS"

while true; do
  case "$1" in
    --genome)
      GENOME_FASTA="$2"
      shift 2
      ;;
    --rnaseq)
      SAMPLE_LIST="$2"
      shift 2
      ;;
    --help)
      usage
      ;;
    --) # End of options
      shift
      break
      ;;
    *)
      echo "Internal error: Unrecognized option '$1'" >&2
      usage
      ;;
  esac
done

# --- Validate Input Arguments ---
if [ -z "${GENOME_FASTA}" ] || [ -z "${SAMPLE_LIST}" ]; then
  echo "Error: Both --genome and --rnaseq arguments are required." >&2
  usage
fi

if [ ! -f "${GENOME_FASTA}" ]; then
  echo "Error: Genome FASTA file not found at ${GENOME_FASTA}" >&2
  exit 1
fi
if [ ! -f "${SAMPLE_LIST}" ]; then
  echo "Error: Sample list file not found at ${SAMPLE_LIST}" >&2
  exit 1
fi

echo "--- Input Files ---"
echo "  Genome FASTA: ${GENOME_FASTA}"
echo "  Sample List: ${SAMPLE_LIST}"
echo ""

# --- 1. Create Output Directories ---
echo "--- Step 1: Creating necessary output directories ---"
mkdir -p "${REF_GENOME_DIR}"
mkdir -p "${BAM_PROCESSING_DIR}"
mkdir -p "${SPLIT_BAM_DIR}"
mkdir -p "${GVCF_DIR}"
mkdir -p "${TEMP_GATK_DIR}"
echo "All output directories created:"
echo "  - Reference Genome: ${REF_GENOME_DIR}"
echo "  - Processed BAMs: ${BAM_PROCESSING_DIR}"
echo "  - Split BAMs: ${SPLIT_BAM_DIR}"
echo "  - GVCFs: ${GVCF_DIR}"
echo "  - GATK Temp: ${TEMP_GATK_DIR}"
echo ""

# --- 2. Prepare Reference Genome for GATK ---
echo "--- Step 2: Preparing Reference Genome for GATK ---"
# GATK requires a .dict file and samtools faidx index (.fai) for the reference.
# These are generated once for the reference genome.
# Copy the genome FASTA to the new REF_GENOME_DIR and work with it there.
GENOME_FASTA_BASENAME=$(basename "${GENOME_FASTA}")
REF_GENOME_PATH="${REF_GENOME_DIR}/${GENOME_FASTA_BASENAME}"

echo "  Copying genome FASTA to ${REF_GENOME_DIR}..."
cp "${GENOME_FASTA}" "${REF_GENOME_PATH}" \
  || { echo "ERROR: Failed to copy genome FASTA. Exiting."; exit 1; }

echo "  Creating sequence dictionary for ${REF_GENOME_PATH}..."
gatk CreateSequenceDictionary -R "${REF_GENOME_PATH}" -O "${REF_GENOME_PATH%.*}.dict" \
  || { echo "ERROR: GATK CreateSequenceDictionary failed. Exiting."; exit 1; }
echo "  Indexing reference genome with samtools faidx..."
samtools faidx "${REF_GENOME_PATH}" \
  || { echo "ERROR: samtools faidx failed. Exiting."; exit 1; }
echo "Reference genome prepared for GATK."
echo ""

# --- 3. Process BAM files for GATK and Call Variants ---
echo "--- Step 3: Processing BAM files for GATK and Calling Variants ---"
echo "  Processing samples from: ${SAMPLE_LIST}"
time while read -r sample_id input_bam_path; do
  # Skip empty lines or lines commented out with '#'
  if [[ -z "$sample_id" || "$sample_id" =~ ^# ]]; then continue; fi

  echo "  Processing sample: ${sample_id}"

  # Define output paths for intermediate files and final GVCF
  UNIQ_BAM="${BAM_PROCESSING_DIR}/${sample_id}.uniq.bam"
  RG_BAM="${BAM_PROCESSING_DIR}/${sample_id}.RG.bam"
  MARKDUP_BAM="${BAM_PROCESSING_DIR}/${sample_id}.markdup.bam"
  MARKDUP_METRICS="${BAM_PROCESSING_DIR}/${sample_id}.markdup_metrics.txt"
  SPLIT_BAM="${SPLIT_BAM_DIR}/${sample_id}.split.bam"
  GVCF_FILE="${GVCF_DIR}/${sample_id}.g.vcf.gz"

  # Ensure the input BAM exists
  if [ ! -f "${input_bam_path}" ]; then
      echo "ERROR: Input BAM not found for ${sample_id} at ${input_bam_path}. Skipping this sample." >&2
      continue
  fi

  # --- 3.1. Index the sorted BAM if it's not already indexed ---
  # Check if .bai file exists, if not, index it.
  if [ ! -f "${input_bam_path}.bai" ]; then
    echo "    3.1. Indexing BAM: ${input_bam_path}"
    samtools index "${input_bam_path}" \
      || { echo "ERROR: samtools index for ${sample_id} failed. Exiting."; exit 1; }
  else
    echo "    3.1. BAM ${input_bam_path} is already indexed."
  fi

  # --- 3.2. Filter uniquely mapped reads and proper pairs ---
  # -F 4 : exclude unmapped reads (0x4)
  # -f 2 : include reads mapped in a proper pair (0x2)
  # -q 255 : include reads with mapping quality 255 (STAR's unique mapping quality)
  echo "    3.2. Filtering uniquely mapped, properly paired reads to: ${UNIQ_BAM}"
  samtools view -@ "${SAMTOOLS_VIEW_THREADS}" -b -h -F 4 -f 2 -q 255 "${input_bam_path}" > "${UNIQ_BAM}" \
    || { echo "ERROR: samtools view (unique filter) for ${sample_id} failed. Exiting."; exit 1; }
  samtools index "${UNIQ_BAM}" \
    || { echo "ERROR: samtools index for unique BAM ${sample_id} failed. Exiting."; exit 1; }

  # --- 3.3. Add or Replace Read Groups ---
  echo "    3.3. Adding/Replacing Read Groups to: ${RG_BAM}"
  gatk AddOrReplaceReadGroups \
    --I "${UNIQ_BAM}" \
    --O "${RG_BAM}" \
    --RGID "${sample_id}" \
    --RGPL ILLUMINA \
    --RGLB "${sample_id}_lib" \
    --RGPU "${sample_id}_pu" \
    --RGSM "${sample_id}" \
    --CREATE_INDEX true \
    || { echo "ERROR: GATK AddOrReplaceReadGroups for ${sample_id} failed. Exiting."; exit 1; }

  # --- 3.4. Mark Duplicates ---
  echo "    3.4. Marking Duplicates to: ${MARKDUP_BAM}"
  gatk --java-options "${MARKDUP_JAVA_OPTIONS}" MarkDuplicates \
    -CREATE_INDEX true \
    -VALIDATION_STRINGENCY SILENT \
    --READ_NAME_REGEX null \
    -I "${RG_BAM}" \
    -M "${MARKDUP_METRICS}" \
    -O "${MARKDUP_BAM}" \
    || { echo "ERROR: GATK MarkDuplicates for ${sample_id} failed. Exiting."; exit 1; }

  # --- 3.5. SplitNCigarReads (for RNA-seq) ---
  # This step is crucial for proper variant calling in RNA-seq data by GATK,
  # as it handles splice junctions.
  echo "    3.5. Splitting N Cigar Reads to: ${SPLIT_BAM}"
  gatk --java-options "${SPLITNCIGAR_JAVA_OPTIONS}" SplitNCigarReads \
    --R "${REF_GENOME_PATH}" \
    --I "${MARKDUP_BAM}" \
    --O "${SPLIT_BAM}" \
    || { echo "ERROR: GATK SplitNCigarReads for ${sample_id} failed. Exiting."; exit 1; }

  # --- 3.6. HaplotypeCaller (GVCF mode) ---
  # Generates a GVCF file per sample. These can then be jointly genotyped
  # across all samples in a later step.
  echo "    3.6. Calling variants with HaplotypeCaller (GVCF) to: ${GVCF_FILE}"
  gatk --java-options "${HAPLOTYPECALLER_JAVA_OPTIONS}" HaplotypeCaller \
    --R "${REF_GENOME_PATH}" \
    --stand-call-conf 20 \
    --dont-use-soft-clipped-bases \
    --I "${SPLIT_BAM}" \
    --O "${GVCF_FILE}" \
    -ERC GVCF \
    || { echo "ERROR: GATK HaplotypeCaller for ${sample_id} failed. Exiting."; exit 1; }

  echo "  GATK processing for ${sample_id} completed."
done < "${SAMPLE_LIST}"
echo "All samples completed GATK processing."
echo ""

echo "--- GATK RNA-seq Variant Calling Pipeline Completed Successfully! ---"
echo "Prepared reference genome files are in: ${REF_GENOME_DIR}"
echo "Processed BAM files are in: ${BAM_PROCESSING_DIR}"
echo "Split N Cigar BAM files are in: ${SPLIT_BAM_DIR}"
echo "GVCF files are in: ${GVCF_DIR}"
echo "You can now proceed with joint genotyping (e.g., GATK GenomicsDBImport and GenotypeGVCFs)."
echo ""
