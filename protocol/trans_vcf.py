import sys
import gzip
import argparse
import os

def process_vcf(input_vcf_path: str, output_vcf_path: str):
    """
    Processes a VCF file (supports .gz compression), modifying its format
    to meet the specified requirements.

    Args:
        input_vcf_path (str): Path to the input VCF file (can be .vcf or .vcf.gz).
        output_vcf_path (str): Path to the output VCF file (can be .vcf or .vcf.gz).
    """
    
    # Determine the appropriate file opening function based on file extension
    _open_infile = gzip.open if input_vcf_path.endswith('.gz') else open
    _open_outfile = gzip.open if output_vcf_path.endswith('.gz') else open

    try:
        # Open files in text mode ('rt' for read text, 'wt' for write text)
        with _open_infile(input_vcf_path, 'rt') as infile, \
             _open_outfile(output_vcf_path, 'wt') as outfile:
            
            for line in infile:
                line = line.strip() # Remove leading/trailing whitespace, including newline characters

                if not line: # Skip empty lines
                    continue

                if line.startswith('##'):
                    # Preserve all meta-information lines
                    outfile.write(line + '\n')
                elif line.startswith('#CHROM'):
                    # Process the header line
                    # Original columns: #CHROM POS ID REF ALT QUAL FILTER INFO FORMAT Sample1 Sample2 ...
                    # Desired columns:  #CHROM POS ID REF ALT QUAL FILTER INFO FORMAT ASE
                    parts = line.split('\t')
                    new_header_parts = parts[:8] # Take the first 8 columns (#CHROM to INFO)
                    new_header_parts.append('FORMAT') # Add the FORMAT column header
                    new_header_parts.append('ASE')    # Add the new ASE column header
                    outfile.write('\t'.join(new_header_parts) + '\n')
                else:
                    # Process data lines
                    # Original data: CHROM POS ID REF ALT QUAL FILTER INFO FORMAT Sample1_data Sample2_data ...
                    # Desired data:  CHROM POS ID REF ALT QUAL FILTER INFO GT 0/1
                    parts = line.split('\t')
                    
                    # Ensure the line has enough columns (at least up to INFO)
                    if len(parts) < 8:
                        sys.stderr.write(f"Warning: Skipping malformed line (fewer than 8 columns): {line}\n")
                        continue

                    # Take the first 8 original data parts (CHROM to INFO)
                    output_parts = parts[:8]
                    output_parts.append('GT')  # New content for the FORMAT column
                    output_parts.append('0/1') # New content for the ASE column

                    outfile.write('\t'.join(output_parts) + '\n')
        
        print(f"VCF file successfully processed. Output saved to: {output_vcf_path}")

    except FileNotFoundError:
        print(f"Error: Input file '{input_vcf_path}' not found. Please check the file path.")
    except Exception as e:
        print(f"An error occurred while processing the file: {e}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Process a VCF file to modify its format. "
                    "It keeps the first 8 columns (CHROM to INFO), "
                    "sets the FORMAT column to 'GT', and adds a new 'ASE' column with '0/1'."
    )
    parser.add_argument(
        '--invcf', 
        type=str, 
        required=True, 
        help='Path to the input VCF file (e.g., your_input.vcf or your_input.vcf.gz).'
    )
    parser.add_argument(
        '--outvcf', 
        type=str, 
        required=True, 
        help='Path to the output VCF file (e.g., your_output.vcf or your_output.vcf.gz).'
    )

    args = parser.parse_args()

    process_vcf(args.invcf, args.outvcf)
