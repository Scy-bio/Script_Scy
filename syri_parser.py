import argparse
import os

def parse_syri_file(input_file, output_format, output_type, output_filename):
    """
    Parsates a syri.out file and outputs selected types in BED or GFF format.

    Args:
        input_file (str): Path to the input syri.out file.
        output_format (str): Desired output format ('bed' or 'gff').
        output_type (str): Type of entries to extract (e.g., 'NOTAL', 'HDR').
        output_filename (str): Name of the output file.
    """
    if output_format.lower() not in ['bed', 'gff']:
        print("Error: --out must be 'bed' or 'gff'.")
        return

    # Prepare a list of valid types to check against
    valid_types = ['SYNAL', 'TRANSAL', 'INVTRAL', 'INVDPAL', 'INVAL', 'DUPAL', 'HDR', 'NOTAL']
    if output_type.upper() not in valid_types:
        print(f"Error: --type must be one of {', '.join(valid_types)}.")
        return

    output_lines = []

    # Add a check for file existence
    if not os.path.exists(input_file):
        print(f"Error: Input file '{input_file}' not found.")
        return

    with open(input_file, 'r') as infile:
        for line in infile:
            parts = line.strip().split('\t')
            # Ensure the line has enough columns
            if len(parts) < 11:
                continue

            # Check the 11th column (index 10) for the specified type
            if parts[10].upper() == output_type.upper():
                # Extract coordinates from column 1-3 and 6-8
                # Column indices are 0-based
                # Column 1 (index 0): HiC_scaffold_1
                # Column 2 (index 1): Start
                # Column 3 (index 2): End

                # Column 6 (index 5): HiC_scaffold_2
                # Column 7 (index 6): Start
                # Column 8 (index 7): End

                # Process first set of coordinates (reference)
                chrom1 = parts[0]
                start1 = None
                end1 = None
                # Only attempt to parse if chrom1 is not '-'
                if chrom1 != '-':
                    try:
                        start1 = int(parts[1])
                        end1 = int(parts[2])
                    except ValueError:
                        # If parsing fails, treat this set of coordinates as invalid
                        print(f"Warning: Could not parse reference coordinates in line: {line.strip()}. Skipping reference output for this line.")
                        chrom1 = '-' # Mark as invalid so it's not outputted

                # Process second set of coordinates (query)
                chrom2 = parts[5]
                start2 = None
                end2 = None
                # Only attempt to parse if chrom2 is not '-'
                if chrom2 != '-':
                    try:
                        start2 = int(parts[6])
                        end2 = int(parts[7])
                    except ValueError:
                        # If parsing fails, treat this set of coordinates as invalid
                        print(f"Warning: Could not parse query coordinates in line: {line.strip()}. Skipping query output for this line.")
                        chrom2 = '-' # Force chrom2 to '-' if parsing fails, so it's not outputted

                # Now, add to output_lines only if valid coordinates exist
                if output_format.lower() == 'bed':
                    if chrom1 != '-':
                        output_lines.append(f"{chrom1}\t{start1 - 1}\t{end1}\t{output_type}\t0\t.")
                    if chrom2 != '-':
                        output_lines.append(f"{chrom2}\t{start2 - 1}\t{end2}\t{output_type}\t0\t.")
                elif output_format.lower() == 'gff':
                    if chrom1 != '-':
                        output_lines.append(f"{chrom1}\tsyri\t{output_type}\t{start1}\t{end1}\t.\t.\t.\tID={output_type}_{chrom1}_{start1}_{end1}")
                    if chrom2 != '-':
                        output_lines.append(f"{chrom2}\tsyri\t{output_type}\t{start2}\t{end2}\t.\t.\t.\tID={output_type}_{chrom2}_{start2}_{end2}")

    # Write to the output file
    if output_filename:
        try:
            with open(output_filename, 'w') as outfile:
                for line in output_lines:
                    outfile.write(line + '\n')
            print(f"Successfully extracted '{output_type}' entries to {output_filename} in {output_format.upper()} format.")
        except IOError as e:
            print(f"Error: Could not write to output file '{output_filename}': {e}")
    else:
        # If no output filename specified, print to stdout (for testing/piping)
        for line in output_lines:
            print(line)

def main():
    parser = argparse.ArgumentParser(description="Parse syri.out file and extract specific entry types into BED or GFF format.")
    parser.add_argument("--input", "-i", dest="input_file", required=True, help="Path to the input syri.out file.")
    parser.add_argument("--out", "-o", dest="output_format", required=True, choices=['bed', 'gff'], help="Desired output format: 'bed' or 'gff'.")
    parser.add_argument("--out_file", "-f", dest="output_filename", help="Name of the output file. If not specified, output will be printed to stdout.")
    parser.add_argument("--type", "-t", dest="output_type", required=True,
                                 choices=['SYNAL', 'TRANSAL', 'INVTRAL', 'INVDPAL', 'INVAL', 'DUPAL', 'HDR', 'NOTAL'],
                                 help="Type of entries to extract (e.g., 'NOTAL', 'HDR', 'SYNAL').")

    args = parser.parse_args()

    parse_syri_file(args.input_file, args.output_format, args.output_type, args.output_filename)

if __name__ == "__main__":
    main()

