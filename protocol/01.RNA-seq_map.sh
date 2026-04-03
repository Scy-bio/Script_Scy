#!/bin/bash
# A comprehensive Bash script for a two-pass STAR transcriptome alignment pipeline.
# Now accepts input file paths as command-line arguments.

# Exit immediately if any command exits with a non-zero status.
set -e

echo "--- Starting STAR Transcriptome Alignment Pipeline ---"
echo "Current working directory: $(pwd)"
echo ""

# --- Configuration Variables ---
# Adjust these values based on your system's available resources (CPU cores, RAM) and data size.

# Number of threads (CPU cores) for STAR genome index generation (Step 2)
GENOME_GENERATE_THREADS=20
# Number of threads for the first pass mapping (Step 4)
PASS1_MAPPING_THREADS=12
# Number of threads for the second pass mapping (Step 6)
PASS2_MAPPING_THREADS=12

# --- Variables for Input Files (will be set by command-line arguments) ---
GENOME_FASTA=""
GTF_FILE=""
SAMPLE_LIST=""

# SJ_OVERHANG parameter for STAR genome index generation.
# This value should typically be (ReadLength - 1).
# For example, if your RNA-seq reads are 150bp long, set this to 149.
# Adjust this value if your reads have a different length.
SJ_OVERHANG=149

# --- Output Directory Structure ---
# All generated output will be stored within this base directory.
STAR_OUTPUT_BASE="test_result/01.RNA-seq_map_results"
INDEX_DIR="${STAR_OUTPUT_BASE}/00.star_index"      # Stores the STAR genome index
PASS1_OUTPUT_DIR="${STAR_OUTPUT_BASE}/01.map_pass1" # Stores first pass alignment logs and SJ.out.tab files
PASS2_OUTPUT_DIR="${STAR_OUTPUT_BASE}/02.map_pass2" # Stores final sorted BAM files and their logs

# --- Usage function ---
usage() {
  echo "Usage: $0 --genome <genome_fasta_file> --gff <gtf_annotation_file> --rnaseq <sample_list_file>"
  echo "Options:"
  echo "  --genome <file>    : Path to the genome FASTA file (e.g., region.fa)"
  echo "  --gff <file>       : Path to the GTF annotation file (e.g., region.gtf)"
  echo "  --rnaseq <file>    : Path to the sample list file (e.g., list)"
  echo "  -h, --help         : Display this help message."
  echo ""
  echo "Example:"
  echo "  $0 --genome data/region.fa --gff data/region.gtf --rnaseq list"
  exit 1
}

# --- 1. Parse Command-line Arguments ---
echo "--- Step 1: Parsing command-line arguments ---"
while [[ "$#" -gt 0 ]]; do
  case "$1" in
    --genome)
      if [ -z "$2" ]; then echo "Error: --genome requires a file path." ; usage; fi
      GENOME_FASTA="$2"
      shift # past argument
      ;;
    --gff)
      if [ -z "$2" ]; then echo "Error: --gff requires a file path." ; usage; fi
      GTF_FILE="$2"
      shift # past argument
      ;;
    --rnaseq)
      if [ -z "$2" ]; then echo "Error: --rnaseq requires a file path." ; usage; fi
      SAMPLE_LIST="$2"
      shift # past argument
      ;;
    -h|--help)
      usage
      ;;
    *)
      echo "Error: Unknown parameter: $1"
      usage
      ;;
  esac
  shift # past argument or value
done

# --- Validate Input Files ---
echo "--- Validating input files ---"
if [ -z "${GENOME_FASTA}" ]; then
  echo "Error: Genome FASTA file (--genome) not specified."
  usage
elif [ ! -f "${GENOME_FASTA}" ]; then
  echo "Error: Genome FASTA file not found at: ${GENOME_FASTA}"
  usage
fi

if [ -z "${GTF_FILE}" ]; then
  echo "Error: GTF annotation file (--gff) not specified."
  usage
elif [ ! -f "${GTF_FILE}" ]; then
  echo "Error: GTF annotation file not found at: ${GTF_FILE}"
  usage
fi

if [ -z "${SAMPLE_LIST}" ]; then
  echo "Error: Sample list file (--rnaseq) not specified."
  usage
elif [ ! -f "${SAMPLE_LIST}" ]; then
  echo "Error: Sample list file not found at: ${SAMPLE_LIST}"
  usage
fi
echo "All required input files found."
echo ""

# --- 2. Create Output Directories ---
echo "--- Step 2: Creating necessary output directories ---"
mkdir -p "${INDEX_DIR}"
mkdir -p "${PASS1_OUTPUT_DIR}"
mkdir -p "${PASS2_OUTPUT_DIR}"
echo "Output directories created:"
echo "  - Index: ${INDEX_DIR}"
echo "  - Pass 1 Maps: ${PASS1_OUTPUT_DIR}"
echo "  - Pass 2 Maps: ${PASS2_OUTPUT_DIR}"
echo ""

# --- 3. Generate STAR Genome Index ---
echo "--- Step 3: Generating STAR Genome Index ---"
echo "  Genome FASTA: ${GENOME_FASTA}"
echo "  GTF Annotation: ${GTF_FILE}"
echo "  Index will be built in: ${INDEX_DIR}"
echo "  This step can take a significant amount of time and memory depending on genome size."
STAR \
  --runMode genomeGenerate \
  --runThreadN "${GENOME_GENERATE_THREADS}" \
  --genomeDir "${INDEX_DIR}" \
  --genomeFastaFiles "${GENOME_FASTA}" \
  --sjdbGTFfile "${GTF_FILE}" \
  --sjdbOverhang "${SJ_OVERHANG}" \
  || { echo "ERROR: STAR genome generation failed. Please check the input files and logs. Exiting."; exit 1; }
echo "STAR genome index generation completed successfully."
echo ""

# --- 4. First Pass STAR Mapping for SJ.out.tab files ---
echo "--- Step 4: Performing First Pass STAR Mapping ---"
echo "  Processing samples from: ${SAMPLE_LIST}"
echo "  Output for this pass will be in: ${PASS1_OUTPUT_DIR}"
time while read -r sample_id r1_file r2_file; do
  # Skip empty lines or lines commented out with '#'
  if [[ -z "$sample_id" || "$sample_id" =~ ^# ]]; then continue; fi

  echo "  Processing sample (Pass 1): ${sample_id}"
  # FASTQ file paths from the list are used directly. Ensure they are correct.
  FULL_R1_PATH="${r1_file}"
  FULL_R2_PATH="${r2_file}"

  # Define log file for this specific sample's first pass
  LOG_FILE="${PASS1_OUTPUT_DIR}/log.STAR_pass1.${sample_id}.txt"

  echo "    Input R1: ${FULL_R1_PATH}, R2: ${FULL_R2_PATH}"
  echo "    Detailed log: ${LOG_FILE}"

  STAR \
    --genomeDir "${INDEX_DIR}" \
    --runThreadN "${PASS1_MAPPING_THREADS}" \
    --readFilesIn "${FULL_R1_PATH}" "${FULL_R2_PATH}" \
    --readFilesCommand zcat \
    --outSAMunmapped Within \
    --outFileNamePrefix "${PASS1_OUTPUT_DIR}/${sample_id}." \
    > "${LOG_FILE}" 2>&1 \
    || { echo "ERROR: STAR Pass 1 for ${sample_id} failed. Check ${LOG_FILE}. Exiting."; exit 1; }
  echo "  First pass mapping for ${sample_id} completed."
done < "${SAMPLE_LIST}"
echo "All samples completed first pass mapping."
echo ""

# --- 5. Filter Splice Junctions from all SJ.out.tab files ---
echo "--- Step 5: Filtering Splice Junctions ---"
SJ_FILTERED_FILE="${PASS1_OUTPUT_DIR}/SJ.filtered.tab"
echo "  Aggregating and filtering all *.SJ.out.tab files from ${PASS1_OUTPUT_DIR}..."
echo "  Output filtered SJ file: ${SJ_FILTERED_FILE}"
cat "${PASS1_OUTPUT_DIR}"/*.SJ.out.tab | \
  awk '($5 > 0 && $7 > 2 && $6==0)' | \
  cut -f1-6 | \
  sort | \
  uniq > "${SJ_FILTERED_FILE}" \
  || { echo "ERROR: Splice junction filtering failed. Exiting."; exit 1; }
echo "Filtered splice junctions saved successfully."
echo ""

# --- 6. Second Pass STAR Mapping using filtered SJs ---
echo "--- Step 6: Performing Second Pass STAR Mapping ---"
echo "  Using filtered splice junctions from: ${SJ_FILTERED_FILE}"
echo "  Final sorted BAM files will be in: ${PASS2_OUTPUT_DIR}"

echo "  Delete the results of First Pass STAR Mapping"
rm "${PASS1_OUTPUT_DIR}"/*.sam

time while read -r sample_id r1_file r2_file; do
  # Skip empty lines or lines commented out with '#'
  if [[ -z "$sample_id" || "$sample_id" =~ ^# ]]; then continue; fi

  echo "  Processing sample (Pass 2): ${sample_id}"
  # FASTQ file paths from the list are used directly. Ensure they are correct.
  FULL_R1_PATH="${r1_file}"
  FULL_R2_PATH="${r2_file}"

  # Define log file for this specific sample's second pass
  LOG_FILE="${PASS2_OUTPUT_DIR}/log.STAR_pass2.${sample_id}.txt"

  echo "    Input R1: ${FULL_R1_PATH}, R2: ${FULL_R2_PATH}"
  echo "    Detailed log: ${LOG_FILE}"

  STAR \
    --genomeDir "${INDEX_DIR}" \
    --runThreadN "${PASS2_MAPPING_THREADS}" \
    --sjdbFileChrStartEnd "${SJ_FILTERED_FILE}" \
    --readFilesIn "${FULL_R1_PATH}" "${FULL_R2_PATH}" \
    --readFilesCommand zcat \
    --outSAMmapqUnique 255 \
    --outSAMtype BAM SortedByCoordinate \
    --outBAMsortingThreadN "${PASS2_MAPPING_THREADS}" \
    --outFileNamePrefix "${PASS2_OUTPUT_DIR}/${sample_id}." \
    > "${LOG_FILE}" 2>&1 \
    || { echo "ERROR: STAR Pass 2 for ${sample_id} failed. Check ${LOG_FILE}. Exiting."; exit 1; }
  echo "  Second pass mapping for ${sample_id} completed."
done < "${SAMPLE_LIST}"
echo "All samples completed second pass mapping."
echo ""

echo "--- STAR Alignment Pipeline Completed Successfully! ---"
echo "Final sorted BAM files are located in: ${PASS2_OUTPUT_DIR}"
echo "You can now proceed with downstream analysis (e.g., counting reads, variant calling)."
echo ""
