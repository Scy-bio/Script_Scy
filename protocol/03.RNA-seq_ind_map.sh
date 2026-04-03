#!/bin/bash
# Script to build an individual-specific genome and perform two-pass RNA-seq alignment with STAR.

# Exit immediately if any command exits with a non-zero status.
set -e

echo "--- Starting Individual Genome Construction and RNA-seq Alignment Pipeline ---"
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

# SJ_OVERHANG parameter for STAR genome index generation.
# This value should typically be (ReadLength - 1).
# For example, if your RNA-seq reads are 150bp long, set this to 149.
# Adjust this value if your reads have a different length.
SJ_OVERHANG=149

# --- Variables for Input Files (will be set by command-line arguments) ---
ORIGINAL_GENOME_FASTA=""  # The reference genome used for variant calling
GTF_FILE=""               # GTF annotation file
SAMPLE_LIST=""            # List of RNA-seq samples (FASTQ files)
SNP2ASE_VCF=""            # T1.snp2ase.vcf.gz path
SNP2MASK_VCF=""           # T1.snp2mask.vcf.gz path
COHORT_ID=""              # The ID for the individual/cohort (e.g., T1)

# --- Output Directory Structure ---
OUTPUT_BASE="test_result/03.ind_map_results"
INDIVIDUAL_GENOME_DIR="${OUTPUT_BASE}/00.individual_genome" # Stores the new individual genome
STAR_INDEX_DIR="${OUTPUT_BASE}/01.star_index"               # Stores the STAR genome index for the individual genome
PASS1_OUTPUT_DIR="${OUTPUT_BASE}/02.map_pass1"              # Stores first pass alignment logs and SJ.out.tab files
PASS2_OUTPUT_DIR="${OUTPUT_BASE}/03.map_pass2"              # Stores final sorted BAM files and their logs

# --- Usage function ---
usage() {
  echo "Usage: $0 --original_genome <ref.fa> --gff <gtf.gtf> --rnaseq <sample_list> --snp2ase <snp2ase.vcf.gz> --snp2mask <snp2mask.vcf.gz> --cohort_id <ID>"
  echo ""
  echo "Options:"
  echo "  --original_genome <file> : Path to the original reference genome FASTA file."
  echo "  --gff <file>             : Path to the GTF annotation file."
  echo "  --rnaseq <file>          : Path to the sample list file (SampleID<tab>R1.fq.gz<tab>R2.fq.gz)."
  echo "  --snp2ase <file>         : Path to the T1.snp2ase.vcf.gz file (used for FastaAlternateReferenceMaker -V)."
  echo "  --snp2mask <file>        : Path to the T1.snp2mask.vcf.gz file (for FastaAlternateReferenceMaker --snp-mask)."
  echo "  --cohort_id <ID>         : Identifier for the individual/cohort (e.g., T1). Used for naming outputs."
  echo "  -h, --help               : Display this help message."
  echo ""
  echo "Example:"
  echo "  $0 --original_genome data/original_ref.fa \\"
  echo "     --gff data/annotation.gtf \\"
  echo "     --rnaseq my_samples.list \\"
  echo "     --snp2ase gatk_joint_calling_results/05.final_vcf/T1.snp2ase.vcf.gz \\"
  echo "     --snp2mask gatk_joint_calling_results/05.final_vcf/T1.snp2mask.vcf.gz \\"
  echo "     --cohort_id T1"
  exit 1
}

# --- Parse Command-line Arguments ---
echo "--- Step 0: Parsing command-line arguments ---"
OPTIONS=$(getopt -o "h" -l original_genome:,gff:,rnaseq:,snp2ase:,snp2mask:,cohort_id:,help -- "$@")

if [ $? -ne 0 ]; then
  echo "Error: Failed to parse command-line arguments." >&2
  usage
fi

eval set -- "$OPTIONS"

while true; do
  case "$1" in
    --original_genome)
      ORIGINAL_GENOME_FASTA="$2"
      shift 2
      ;;
    --gff)
      GTF_FILE="$2"
      shift 2
      ;;
    --rnaseq)
      SAMPLE_LIST="$2"
      shift 2
      ;;
    --snp2ase)
      SNP2ASE_VCF="$2"
      shift 2
      ;;
    --snp2mask)
      SNP2MASK_VCF="$2"
      shift 2
      ;;
    --cohort_id)
      COHORT_ID="$2"
      shift 2
      ;;
    -h|--help)
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

# --- Validate Input Files ---
echo "--- Validating input files and arguments ---"
if [ -z "${ORIGINAL_GENOME_FASTA}" ] || [ -z "${GTF_FILE}" ] || [ -z "${SAMPLE_LIST}" ] || \
   [ -z "${SNP2ASE_VCF}" ] || [ -z "${SNP2MASK_VCF}" ] || [ -z "${COHORT_ID}" ]; then
  echo "Error: All required arguments (--original_genome, --gff, --rnaseq, --snp2ase, --snp2mask, --cohort_id) must be specified."
  usage
fi

for file in "${ORIGINAL_GENOME_FASTA}" "${GTF_FILE}" "${SAMPLE_LIST}" "${SNP2ASE_VCF}" "${SNP2MASK_VCF}"; do
  if [ ! -f "${file}" ]; then
    echo "Error: Input file not found: ${file}" >&2
    exit 1
  fi
done
echo "All required input files and arguments are valid."
echo ""

# --- 1. Create Output Directories ---
echo "--- Step 1: Creating necessary output directories ---"
mkdir -p "${INDIVIDUAL_GENOME_DIR}"
mkdir -p "${STAR_INDEX_DIR}"
mkdir -p "${PASS1_OUTPUT_DIR}"
mkdir -p "${PASS2_OUTPUT_DIR}"
echo "Output directories created:"
echo "  - Individual Genome: ${INDIVIDUAL_GENOME_DIR}"
echo "  - STAR Index: ${STAR_INDEX_DIR}"
echo "  - Pass 1 Maps: ${PASS1_OUTPUT_DIR}"
echo "  - Pass 2 Maps: ${PASS2_OUTPUT_DIR}"
echo ""

# Define the path for the newly generated individual genome
INDIVIDUAL_GENOME_FASTA="${INDIVIDUAL_GENOME_DIR}/${COHORT_ID}.genome.fa"

# --- 2. Construct Individual-Specific Reference Genome ---
echo "--- Step 2: Constructing Individual-Specific Reference Genome using GATK FastaAlternateReferenceMaker ---"
echo "  Original Reference: ${ORIGINAL_GENOME_FASTA}"
echo "  Variants (Alt): ${SNP2ASE_VCF}"
echo "  SNP Mask: ${SNP2MASK_VCF}"
echo "  Output Individual Genome: ${INDIVIDUAL_GENOME_FASTA}"
echo "  NOTE: This step will incorporate variants from ${SNP2ASE_VCF} and mask sites from ${SNP2MASK_VCF}."
echo "  The --snp-mask-priority flag ensures masked sites take precedence."

gatk FastaAlternateReferenceMaker \
  -R "${ORIGINAL_GENOME_FASTA}" \
  --line-width 80 \
  -V "${SNP2MASK_VCF}"  \
  --snp-mask  "${SNP2ASE_VCF}" \
  --snp-mask-priority \
  -O "${INDIVIDUAL_GENOME_FASTA}" \
  || { echo "ERROR: GATK FastaAlternateReferenceMaker failed. Exiting."; exit 1; }
echo "Individual genome construction completed."
echo ""

# --- 3. Correct FASTA Headers (if needed) ---
# This awk command is often necessary because FastaAlternateReferenceMaker can add
# extra text to headers that STAR might not like.
echo "--- Step 3: Correcting FASTA Headers in ${INDIVIDUAL_GENOME_FASTA} ---"
awk -F "[ :]" '{if($1~">")print ">"$2 ; else print$0}' "${INDIVIDUAL_GENOME_FASTA}" > "${INDIVIDUAL_GENOME_FASTA}.tmp" && \
  mv "${INDIVIDUAL_GENOME_FASTA}.tmp" "${INDIVIDUAL_GENOME_FASTA}" \
  || { echo "ERROR: FASTA header correction failed. Exiting."; exit 1; }
echo "FASTA headers corrected."
echo ""

# --- 4. Prepare Individual Genome for GATK/Samtools ---
echo "--- Step 4: Preparing Individual Genome for GATK/Samtools ---"
# GATK/Samtools require a .dict file and samtools faidx index (.fai) for the reference.
echo "  Creating sequence dictionary for ${INDIVIDUAL_GENOME_FASTA}..."
rm "${INDIVIDUAL_GENOME_FASTA%.*}.dict" "${INDIVIDUAL_GENOME_FASTA%.*}.fa.fai" 
gatk CreateSequenceDictionary -R "${INDIVIDUAL_GENOME_FASTA}" -O "${INDIVIDUAL_GENOME_FASTA%.*}.dict" \
  || { echo "ERROR: GATK CreateSequenceDictionary for individual genome failed. Exiting."; exit 1; }
echo "  Indexing individual genome with samtools faidx..."
samtools faidx "${INDIVIDUAL_GENOME_FASTA}" \
  || { echo "ERROR: samtools faidx for individual genome failed. Exiting."; exit 1; }
echo "Individual genome prepared for GATK/Samtools."
echo ""

# --- 5. Generate STAR Genome Index for Individual Genome ---
echo "--- Step 5: Generating STAR Genome Index for ${COHORT_ID} ---"
echo "  Genome FASTA: ${INDIVIDUAL_GENOME_FASTA}"
echo "  GTF Annotation: ${GTF_FILE}"
echo "  Index will be built in: ${STAR_INDEX_DIR}"
echo "  This step can take a significant amount of time and memory depending on genome size."
STAR \
  --runMode genomeGenerate \
  --runThreadN "${GENOME_GENERATE_THREADS}" \
  --genomeDir "${STAR_INDEX_DIR}" \
  --genomeFastaFiles "${INDIVIDUAL_GENOME_FASTA}" \
  --sjdbGTFfile "${GTF_FILE}" \
  --sjdbOverhang "${SJ_OVERHANG}" \
  || { echo "ERROR: STAR genome generation for individual genome failed. Please check the input files and logs. Exiting."; exit 1; }
echo "STAR genome index generation for individual genome completed successfully."
echo ""

# --- 6. First Pass STAR Mapping for SJ.out.tab files ---
echo "--- Step 6: Performing First Pass STAR Mapping to Individual Genome ---"
echo "  Processing samples from: ${SAMPLE_LIST}"
echo "  Output for this pass will be in: ${PASS1_OUTPUT_DIR}"
time while read -r sample_id r1_file r2_file; do
  # Skip empty lines or lines commented out with '#'
  if [[ -z "$sample_id" || "$sample_id" =~ ^# ]]; then continue; fi

  # Validate FASTQ files
  if [ ! -f "${r1_file}" ]; then
    echo "ERROR: R1 FASTQ file not found for ${sample_id}: ${r1_file}. Skipping." >&2
    continue
  fi
  if [ ! -f "${r2_file}" ]; then
    echo "ERROR: R2 FASTQ file not found for ${sample_id}: ${r2_file}. Skipping." >&2
    continue
  fi

  echo "  Processing sample (Pass 1): ${sample_id}"
  # Define log file for this specific sample's first pass
  LOG_FILE="${PASS1_OUTPUT_DIR}/log.STAR_pass1.${sample_id}.txt"

  echo "    Input R1: ${r1_file}, R2: ${r2_file}"
  echo "    Detailed log: ${LOG_FILE}"

  STAR \
    --genomeDir "${STAR_INDEX_DIR}" \
    --runThreadN "${PASS1_MAPPING_THREADS}" \
    --readFilesIn "${r1_file}" "${r2_file}" \
    --readFilesCommand zcat \
    --outSAMunmapped Within \
    --outFileNamePrefix "${PASS1_OUTPUT_DIR}/${sample_id}." \
    > "${LOG_FILE}" 2>&1 \
    || { echo "ERROR: STAR Pass 1 for ${sample_id} failed. Check ${LOG_FILE}. Exiting."; exit 1; }
  echo "  First pass mapping for ${sample_id} completed."
done < "${SAMPLE_LIST}"
echo "All samples completed first pass mapping."
echo ""

# --- 7. Filter Splice Junctions from all SJ.out.tab files ---
echo "--- Step 7: Filtering Splice Junctions ---"
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

# --- 8. Second Pass STAR Mapping using filtered SJs ---
echo "--- Step 8: Performing Second Pass STAR Mapping to Individual Genome ---"
echo "  Using filtered splice junctions from: ${SJ_FILTERED_FILE}"
echo "  Final sorted BAM files will be in: ${PASS2_OUTPUT_DIR}"

echo "  Deleting SAM files from First Pass STAR Mapping (to save space)..."
# Using find with -delete is safer than rm *.sam
find "${PASS1_OUTPUT_DIR}" -name "*.sam" -delete
echo "  First pass SAM files deleted."
echo ""

time while read -r sample_id r1_file r2_file; do
  # Skip empty lines or lines commented out with '#'
  if [[ -z "$sample_id" || "$sample_id" =~ ^# ]]; then continue; fi

  # Validate FASTQ files again (redundant but safe)
  if [ ! -f "${r1_file}" ]; then
    echo "ERROR: R1 FASTQ file not found for ${sample_id}: ${r1_file}. Skipping." >&2
    continue
  fi
  if [ ! -f "${r2_file}" ]; then
    echo "ERROR: R2 FASTQ file not found for ${sample_id}: ${r2_file}. Skipping." >&2
    continue
  fi

  echo "  Processing sample (Pass 2): ${sample_id}"
  # Define log file for this specific sample's second pass
  LOG_FILE="${PASS2_OUTPUT_DIR}/log.STAR_pass2.${sample_id}.txt"

  echo "    Input R1: ${r1_file}, R2: ${r2_file}"
  echo "    Detailed log: ${LOG_FILE}"

  STAR \
    --genomeDir "${STAR_INDEX_DIR}" \
    --runThreadN "${PASS2_MAPPING_THREADS}" \
    --sjdbFileChrStartEnd "${SJ_FILTERED_FILE}" \
    --readFilesIn "${r1_file}" "${r2_file}" \
    --readFilesCommand zcat \
    --outSAMmapqUnique 255 \
    --outSAMtype BAM SortedByCoordinate \
    --outBAMsortingThreadN "${PASS2_MAPPING_THREADS}" \
    --outFileNamePrefix "${PASS2_OUTPUT_DIR}/${sample_id}." \
    > "${LOG_FILE}" 2>&1 \
    || { echo "ERROR: STAR Pass 2 for ${sample_id} failed. Check ${LOG_FILE}. Exiting."; exit 1; }

  # Index the final BAM file
  echo "  Indexing final BAM file for ${sample_id}..."
  samtools index "${PASS2_OUTPUT_DIR}/${sample_id}.Aligned.sortedByCoord.out.bam" \
    || { echo "WARNING: Samtools indexing for ${sample_id}.Aligned.sortedByCoord.out.bam failed. Continuing."; }

  echo "  Second pass mapping for ${sample_id} completed."
done < "${SAMPLE_LIST}"
echo "All samples completed second pass mapping."
echo ""

echo "--- Individual Genome Construction and RNA-seq Alignment Pipeline Completed Successfully! ---"
echo "Individual genome FASTA: ${INDIVIDUAL_GENOME_FASTA}"
echo "Final sorted BAM files are located in: ${PASS2_OUTPUT_DIR}"
echo "You can now proceed with downstream analysis (e.g., counting reads, variant calling on the new genome)."
echo ""
