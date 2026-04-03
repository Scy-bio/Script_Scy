#!/bin/bash
# A script for GATK joint genotyping and SNP filtration.
# Takes GVCF files from individual samples and a reference genome as input.
# Now includes options to specify the path to the filter_snps.py script,
# and DP and GP thresholds for the Python script.

# Exit immediately if any command exits with a non-zero status.
set -e

echo "--- Starting GATK Joint Genotyping and SNP Filtration Pipeline ---"
echo "Current working directory: $(pwd)"
echo ""

# --- Default Configuration Variables ---
# These can be overridden by command-line arguments or adjusted manually.

# Input file paths (will be set by command-line arguments)
GENOME_FASTA=""          # Path to the reference genome FASTA file
GVCF_LIST_FILE=""        # Path to the list file containing GVCF paths
FILTER_SNPS_PYTHON_SCRIPT="" # Path to the custom Python script filter_snps.py

# Default thresholds for filter_snps.py
FILTER_DP_THRESHOLD=10   # Default minimum depth
FILTER_GP_THRESHOLD=0.5  # Default minimum genotype quality

# GATK Java Options (adjust -Xmx for memory and -XX:ParallelGCThreads for cores)
# Ensure -Djava.io.tmpdir points to a directory with sufficient space for temporary files.
# These values are examples; adjust them based on your system's resources.
GATK_JAVA_OPTIONS_COMBINE="-Xmx12g" # For CombineGVCFs
GATK_JAVA_OPTIONS_GENOTYPE="-Xmx12g" # For GenotypeGVCFs
GATK_JAVA_OPTIONS_DEFAULT="-Xmx4g"  # For other GATK tools if specific options aren't set

# --- Output Directory Structure ---
OUTPUT_BASE="test_result/02.SNP_filter_result"
GVCF_COMBINED_DIR="${OUTPUT_BASE}/01.combined_gvcf"
GENOTYPED_VCF_DIR="${OUTPUT_BASE}/02.genotyped_vcf"
FILTERED_VCF_DIR="${OUTPUT_BASE}/03.filtered_vcf"
QC_OUTPUT_DIR="${OUTPUT_BASE}/04.qc_outputs"
FINAL_VCF_DIR="${OUTPUT_BASE}/05.final_vcf"
TEMP_DIR="${OUTPUT_BASE}/tmp_gatk_joint" # Temporary directory for GATK operations

# --- Usage Function ---
usage() {
  echo "Usage: $0 --genome <path/to/ref.fa> --gvcf_list <path/to/list_file> --script <path/to/filter_snps.py> [--dp <min_depth>] [--gp <min_gp>]"
  echo ""
  echo "Options:"
  echo "  --genome      Path to the reference genome FASTA file (e.g., data/ref.fa)."
  echo "  --gvcf_list   Path to the list file of GVCFs. Each line should be:"
  echo "                CohortID<tab>Path/To/Sample.g.vcf.gz"
  echo "  --script      Path to your custom Python script filter_snps.py."
  echo "  --dp          Minimum Depth (DP) threshold (Note: >=20% tissues) for filter_snps.py (default: ${FILTER_DP_THRESHOLD})."
  echo "  --gp          Minimum Genotype Quality (GP) threshold for filter_snps.py (default: ${FILTER_GP_THRESHOLD})."
  echo ""
  echo "Example:"
  echo "  ./run_joint_genotyping.sh --genome data/region.fa --gvcf_list my_cohort.list --script ~/my_scripts/filter_snps.py --dp 15 --gp 0.7"
  exit 1
}

# --- Parse Command Line Arguments ---
OPTIONS=$(getopt -o "" -l genome:,gvcf_list:,script:,dp:,gp: -- "$@")

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
    --gvcf_list)
      GVCF_LIST_FILE="$2"
      shift 2
      ;;
    --script)
      FILTER_SNPS_PYTHON_SCRIPT="$2"
      shift 2
      ;;
    --dp)
      FILTER_DP_THRESHOLD="$2"
      shift 2
      ;;
    --gp)
      FILTER_GP_THRESHOLD="$2"
      shift 2
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
if [ -z "${GENOME_FASTA}" ] || [ -z "${GVCF_LIST_FILE}" ] || [ -z "${FILTER_SNPS_PYTHON_SCRIPT}" ]; then
  echo "Error: --genome, --gvcf_list, and --script arguments are required." >&2
  usage
fi

if [ ! -f "${GENOME_FASTA}" ]; then
  echo "Error: Reference genome FASTA file not found at ${GENOME_FASTA}" >&2
  exit 1
fi
if [ ! -f "${GVCF_LIST_FILE}" ]; then
  echo "Error: GVCF list file not found at ${GVCF_LIST_FILE}" >&2
  exit 1
fi
if [ ! -f "${FILTER_SNPS_PYTHON_SCRIPT}" ]; then
  echo "ERROR: Custom Python script '${FILTER_SNPS_PYTHON_SCRIPT}' not found at the specified path." >&2
  echo "Please ensure the path is correct and the script exists. Exiting." >&2
  exit 1
fi

# Basic validation for DP and GP thresholds (check if numeric)
if ! [[ "${FILTER_DP_THRESHOLD}" =~ ^[0-9]+$ ]]; then
  echo "ERROR: Invalid value for --dp. Please provide an integer (e.g., 10)." >&2
  usage
fi
# Using a regex that allows integers or decimals for GP
if ! [[ "${FILTER_GP_THRESHOLD}" =~ ^[0-9]+(\.[0-9]+)?$ ]]; then
  echo "ERROR: Invalid value for --gp. Please provide a numeric value (e.g., 0.5)." >&2
  usage
fi


echo "--- Input Files ---"
echo "  Reference Genome FASTA: ${GENOME_FASTA}"
echo "  GVCF List File: ${GVCF_LIST_FILE}"
echo "  Python Filter Script: ${FILTER_SNPS_PYTHON_SCRIPT}"
echo "  Filter DP Threshold: ${FILTER_DP_THRESHOLD}"
echo "  Filter GP Threshold: ${FILTER_GP_THRESHOLD}"
echo ""

# --- 1. Create Output Directories ---
echo "--- Step 1: Creating necessary output directories ---"
mkdir -p "${GVCF_COMBINED_DIR}"
mkdir -p "${GENOTYPED_VCF_DIR}"
mkdir -p "${FILTERED_VCF_DIR}"
mkdir -p "${QC_OUTPUT_DIR}"
mkdir -p "${FINAL_VCF_DIR}"
mkdir -p "${TEMP_DIR}"
echo "All output directories created."
echo ""

# --- 2. Prepare Reference Genome for GATK (if not already done) ---
# GATK requires a .dict file and samtools faidx index (.fai) for the reference.
# Check if these exist; if not, create them.
if [ ! -f "${GENOME_FASTA%.*}.dict" ]; then
  echo "--- Step 2: Creating sequence dictionary for ${GENOME_FASTA} ---"
  gatk CreateSequenceDictionary -R "${GENOME_FASTA}" -O "${GENOME_FASTA%.*}.dict" \
    || { echo "ERROR: GATK CreateSequenceDictionary failed. Exiting."; exit 1; }
fi
if [ ! -f "${GENOME_FASTA}.fai" ]; then
  echo "--- Step 2: Indexing reference genome with samtools faidx ---"
  samtools faidx "${GENOME_FASTA}" \
    || { echo "ERROR: samtools faidx failed. Exiting."; exit 1; }
fi
echo "Reference genome prepared for GATK."
echo ""

# --- Extract Cohort ID and GVCF paths ---
# The first column of the list file is assumed to be the Cohort ID
# All GVCFs in the list will be combined for this cohort.
COHORT_ID=$(awk 'NR==1 {print $1}' "${GVCF_LIST_FILE}")
if [ -z "$COHORT_ID" ]; then
    echo "ERROR: Could not extract Cohort ID from ${GVCF_LIST_FILE}. Is the file empty or malformed?" >&2
    exit 1
fi
echo "Processing cohort: ${COHORT_ID}"

# Create the GVCFs list file for CombineGVCFs
GVCF_INPUT_LIST="${TEMP_DIR}/${COHORT_ID}.gvcf.list"
echo "Extracting GVCF paths to ${GVCF_INPUT_LIST}..."
awk '{print $2}' "${GVCF_LIST_FILE}" > "${GVCF_INPUT_LIST}"

# Check if the extracted GVCFs actually exist
while read -r -a line_array; do
    gvcf_file="${line_array[1]}" # The second element is the GVCF path
    if [[ ! "$gvcf_file" =~ \.g\.vcf\.gz$ ]]; then
        echo "WARNING: GVCF file path does not end with .g.vcf.gz: ${gvcf_file}"
    fi
    if [ ! -f "$gvcf_file" ]; then
        echo "ERROR: GVCF file not found: ${gvcf_file}. Please check your list file. Exiting." >&2
        exit 1
    fi
done < "${GVCF_LIST_FILE}"
echo "All GVCFs in the list appear to be valid."
echo ""

# --- 3. Combine GVCFs ---
echo "--- Step 3: Combining GVCFs for Cohort ${COHORT_ID} ---"
COMBINED_GVCF="${GVCF_COMBINED_DIR}/${COHORT_ID}.g.vcf.gz"
echo "  Output combined GVCF: ${COMBINED_GVCF}"

gatk --java-options "${GATK_JAVA_OPTIONS_COMBINE}" CombineGVCFs \
  -R "${GENOME_FASTA}" \
  -V "${GVCF_INPUT_LIST}" \
  -O "${COMBINED_GVCF}" \
  || { echo "ERROR: GATK CombineGVCFs failed for ${COHORT_ID}. Exiting."; exit 1; }
echo "Combined GVCFs successfully."
echo ""

# --- 4. Genotype GVCFs ---
echo "--- Step 4: Genotyping Combined GVCF for Cohort ${COHORT_ID} ---"
GENOTYPED_VCF="${GENOTYPED_VCF_DIR}/${COHORT_ID}.vcf.gz"
echo "  Output genotyped VCF: ${GENOTYPED_VCF}"
gatk --java-options "${GATK_JAVA_OPTIONS_GENOTYPE}" GenotypeGVCFs \
  -R "${GENOME_FASTA}" \
  --variant "${COMBINED_GVCF}" \
  -O "${GENOTYPED_VCF}" \
  || { echo "ERROR: GATK GenotypeGVCFs failed for ${COHORT_ID}. Exiting."; exit 1; }
echo "Genotyping completed successfully."
echo ""

# --- 5. Select SNPs ---
echo "--- Step 5: Selecting Biallelic SNPs for Cohort ${COHORT_ID} ---"
SELECTED_SNP_VCF="${FILTERED_VCF_DIR}/${COHORT_ID}.snp.vcf.gz"
echo "  Output SNP VCF: ${SELECTED_SNP_VCF}"
gatk --java-options "${GATK_JAVA_OPTIONS_DEFAULT}" SelectVariants \
  -select-type SNP \
  --restrict-alleles-to BIALLELIC \
  -V "${GENOTYPED_VCF}" \
  -O "${SELECTED_SNP_VCF}" \
  || { echo "ERROR: GATK SelectVariants (SNP) failed for ${COHORT_ID}. Exiting."; exit 1; }
echo "SNP selection completed."
echo ""

# --- 6. Variant Filtration ---
echo "--- Step 6: Filtering SNPs for Cohort ${COHORT_ID} (FS and QD) ---"
FILTERED_SNP_VCF="${FILTERED_VCF_DIR}/${COHORT_ID}.snps_fld.vcf.gz"
echo "  Output filtered SNP VCF: ${FILTERED_SNP_VCF}"
gatk --java-options "${GATK_JAVA_OPTIONS_DEFAULT}" VariantFiltration \
  --filter-name "FS" \
  --filter "FS > 30.0" \
  --filter-name "QD" \
  --filter "QD < 2.0" \
  -V "${SELECTED_SNP_VCF}" \
  -O "${FILTERED_SNP_VCF}" \
  || { echo "ERROR: GATK VariantFiltration failed for ${COHORT_ID}. Exiting."; exit 1; }
echo "SNP filtration completed."
echo ""

# --- 7. Select Non-Filtered SNPs ---
echo "--- Step 7: Selecting Non-Filtered SNPs for Cohort ${COHORT_ID} ---"
FINAL_FILTERED_SNP_VCF="${FINAL_VCF_DIR}/${COHORT_ID}.snps_fld_sel.vcf.gz"
echo "  Output final selected SNP VCF: ${FINAL_FILTERED_SNP_VCF}"
gatk --java-options "${GATK_JAVA_OPTIONS_DEFAULT}" SelectVariants \
  --exclude-filtered TRUE \
  -V "${FILTERED_VCF_DIR}/${COHORT_ID}.snps_fld.vcf.gz" \
  -O "${FINAL_FILTERED_SNP_VCF}" \
  || { echo "ERROR: GATK SelectVariants (exclude filtered) failed for ${COHORT_ID}. Exiting."; exit 1; }
echo "Non-filtered SNP selection completed."
echo ""

# --- 8. Extract Data for Python Script ---
echo "--- Step 8: Extracting SNP data for Python processing for Cohort ${COHORT_ID} ---"
BCFTOOLS_OUTPUT_TXT="${QC_OUTPUT_DIR}/${COHORT_ID}.output.txt"
echo "  Output BCFtools query result: ${BCFTOOLS_OUTPUT_TXT}"
# Check if bcftools exists and is executable
if ! command -v bcftools &> /dev/null; then
    echo "ERROR: bcftools command not found. Please install bcftools and ensure it's in your PATH. Exiting." >&2
    exit 1
fi
bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT[\t%GT\t%DP]\n' "${FINAL_FILTERED_SNP_VCF}" > "${BCFTOOLS_OUTPUT_TXT}" \
  || { echo "ERROR: bcftools query failed for ${COHORT_ID}. Exiting."; exit 1; }
echo "BCFtools query completed."
echo ""

# --- 9. Run Custom Python Script ---
echo "--- Step 9: Running custom Python script (${FILTER_SNPS_PYTHON_SCRIPT}) for Cohort ${COHORT_ID} ---"
PYTHON_FILTERED_TXT="${QC_OUTPUT_DIR}/${COHORT_ID}.filtered.txt"
PYTHON_SITES_LIST="${QC_OUTPUT_DIR}/${COHORT_ID}.sites.list"
echo "  Output filtered data: ${PYTHON_FILTERED_TXT}"
echo "  Output sites list: ${PYTHON_SITES_LIST}"
# Ensure the python script is executable
if [ ! -x "${FILTER_SNPS_PYTHON_SCRIPT}" ]; then
    echo "WARNING: Python script '${FILTER_SNPS_PYTHON_SCRIPT}' is not executable. Attempting to run with 'python3'." >&2
    python3 "${FILTER_SNPS_PYTHON_SCRIPT}" \
      --in "${BCFTOOLS_OUTPUT_TXT}" \
      --out "${PYTHON_FILTERED_TXT}" \
      --DP "${FILTER_DP_THRESHOLD}" \
      --GP "${FILTER_GP_THRESHOLD}" \
      --list "${PYTHON_SITES_LIST}" \
      || { echo "ERROR: Custom Python script '${FILTER_SNPS_PYTHON_SCRIPT}' failed for ${COHORT_ID}. Exiting."; exit 1; }
else
    "${FILTER_SNPS_PYTHON_SCRIPT}" \
      --in "${BCFTOOLS_OUTPUT_TXT}" \
      --out "${PYTHON_FILTERED_TXT}" \
      --DP "${FILTER_DP_THRESHOLD}" \
      --GP "${FILTER_GP_THRESHOLD}" \
      --list "${PYTHON_SITES_LIST}" \
      || { echo "ERROR: Custom Python script '${FILTER_SNPS_PYTHON_SCRIPT}' failed for ${COHORT_ID}. Exiting."; exit 1; }
fi
echo "Custom Python script completed."
echo ""

# --- 10. Final Variant Selection based on Python Script Output ---
echo "--- Step 10: Final Variant Selection based on Python script output for Cohort ${COHORT_ID} ---"
SNP2MASK_VCF="${FINAL_VCF_DIR}/${COHORT_ID}.snp2mask.vcf.gz"
SNP2ASE_VCF="${FINAL_VCF_DIR}/${COHORT_ID}.snp2ase.vcf.gz"
echo "  Output VCF for masking: ${SNP2MASK_VCF}"
echo "  Output VCF for allele specific expression: ${SNP2ASE_VCF}"

# Select variants for masking
gatk --java-options "${GATK_JAVA_OPTIONS_DEFAULT}" SelectVariants \
  -L "${PYTHON_SITES_LIST}" \
  -V "${FINAL_FILTERED_SNP_VCF}" \
  -O "${SNP2MASK_VCF}" \
  || { echo "ERROR: GATK SelectVariants (snp2mask) failed for ${COHORT_ID}. Exiting."; exit 1; }

# Select variants for allele-specific expression (AF<1.00 implies not monomorphic)
gatk --java-options "${GATK_JAVA_OPTIONS_DEFAULT}" SelectVariants \
  -L "${PYTHON_SITES_LIST}" \
  -select "AF<1.00" \
  -V "${FINAL_FILTERED_SNP_VCF}" \
  -O "${SNP2ASE_VCF}" \
  || { echo "ERROR: GATK SelectVariants (snp2ase) failed for ${COHORT_ID}. Exiting."; exit 1; }
echo "Final variant selection completed."
echo ""

echo "--- GATK Joint Genotyping and SNP Filtration Pipeline Completed Successfully! ---"
echo "All outputs are in the '${OUTPUT_BASE}' directory."
echo "You can now use '${FINAL_VCF_DIR}/${COHORT_ID}.snp2mask.vcf.gz' for masking or"
echo "'${FINAL_VCF_DIR}/${COHORT_ID}.snp2ase.vcf.gz' for allele-specific expression analysis."
echo ""
