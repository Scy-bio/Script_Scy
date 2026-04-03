#!/usr/bin/env python3

import argparse
import subprocess
import os

def run_command(cmd_list, shell=False):
    """Helper function to run shell commands."""
    try:
        if shell:
            # When shell=True, cmd_list should be a single string
            subprocess.run(cmd_list, check=True, text=True, shell=True, capture_output=True)
        else:
            # When shell=False (default), cmd_list should be a list of arguments
            subprocess.run(cmd_list, check=True, text=True, capture_output=True)
    except subprocess.CalledProcessError as e:
        print(f"Error running command: {' '.join(cmd_list) if isinstance(cmd_list, list) else cmd_list}")
        print(f"STDOUT: {e.stdout}")
        print(f"STDERR: {e.stderr}")
        exit(1)

def parse_lengths(filepath):
    """Parses a length file (e.g., A.length) into a dictionary."""
    lengths = {}
    with open(filepath, 'r') as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) == 2:
                try:
                    lengths[parts[0]] = int(parts[1])
                except ValueError:
                    print(f"Warning: Could not parse length for {parts[0]} in {filepath}. Skipping.")
    return lengths

def main():
    parser = argparse.ArgumentParser(
        description="Extracts gene lengths and filters allele pairs based on length ratio."
    )
    parser.add_argument(
        "--hapA",
        required=True,
        help="Path to HapA gene FASTA file (e.g., CskiameaA.gene.fa)"
    )
    parser.add_argument(
        "--hapB",
        required=True,
        help="Path to HapB gene FASTA file (e.g., CskiameaB.gene.fa)"
    )
    parser.add_argument(
        "--allele_table",
        required=True,
        help="Path to the allele table file (e.g., allele.table)"
    )
    parser.add_argument(
        "--filter_length",
        type=float,
        default=1.25,
        help="Filtering threshold for length difference. Keep pairs where (longer_length / shorter_length) < filter_length. Default is 1.25 (25%% difference)."
    )
    parser.add_argument(
        "--out",
        required=True,
        help="Output file name for the filtered allele length table (e.g., allele.length.table)"
    )

    args = parser.parse_args()

    # Step 1: Extract gene length information using seqkit
    print("Step 1: Extracting gene lengths...")
    a_length_file = "A.length"
    b_length_file = "B.length"

    # CORRECTED PART: Use shell=True for redirection or directly redirect output
    # Option 1: Use shell=True (simpler for direct shell commands)
    run_command(f"seqkit fx2tab -l -n {args.hapA} > {a_length_file}", shell=True)
    run_command(f"seqkit fx2tab -l -n {args.hapB} > {b_length_file}", shell=True)

    # Option 2 (more robust Pythonic way without shell=True):
    # with open(a_length_file, 'w') as outfile_a:
    #     subprocess.run(["seqkit", "fx2tab", "-l", "-n", args.hapA], check=True, text=True, stdout=outfile_a)
    # with open(b_length_file, 'w') as outfile_b:
    #     subprocess.run(["seqkit", "fx2tab", "-l", "-n", args.hapB], check=True, text=True, stdout=outfile_b)


    print(f"Lengths extracted to {a_length_file} and {b_length_file}.")

    # Load lengths into dictionaries
    print("Loading gene lengths into memory...")
    a_lengths = parse_lengths(a_length_file)
    b_lengths = parse_lengths(b_length_file)
    print(f"Loaded {len(a_lengths)} lengths for HapA and {len(b_lengths)} for HapB.")

    # Step 2: Merge and filter allele lengths
    print(f"Step 2: Merging and filtering allele lengths (filter_length={args.filter_length})...")
    
    filtered_count = 0
    with open(args.allele_table, 'r') as infile, open(args.out, 'w') as outfile:
        # Write header to output file
        outfile.write("GeneA\tGeneB\tLengthA\tLengthB\n")
        
        for line in infile:
            parts = line.strip().split()
            # The allele table usually has two columns. If 04.result has more, parts[0] and parts[3] might be needed.
            # Based on your example `CsiA_G00000000007       CsiB_G00000000006`, it's just two columns.
            if len(parts) < 2:
                continue # Skip malformed lines or empty lines

            gene_a = parts[0]
            gene_b = parts[1]

            len_a = a_lengths.get(gene_a)
            len_b = b_lengths.get(gene_b)

            if len_a is None:
                print(f"Warning: Length for {gene_a} (from allele table) not found in {a_length_file}. Skipping pair.")
                continue
            if len_b is None:
                print(f"Warning: Length for {gene_b} (from allele table) not found in {b_length_file}. Skipping pair.")
                continue

            shorter_len = min(len_a, len_b)
            longer_len = max(len_a, len_b)

            if shorter_len == 0:
                continue 

            # Filtering condition: longer_length / shorter_length < filter_length
            if (longer_len / shorter_len) < args.filter_length:
                outfile.write(f"{gene_a}\t{gene_b}\t{len_a}\t{len_b}\n")
                filtered_count += 1
    
    print(f"Filtering complete. {filtered_count} allele pairs written to {args.out}.")

    # Clean up temporary length files
    print("Cleaning up temporary files...")
    os.remove(a_length_file)
    os.remove(b_length_file)
    print("Temporary files removed.")

if __name__ == "__main__":
    main()
