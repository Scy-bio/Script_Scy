#!/usr/bin/env python3

import argparse
import subprocess
import os
import sys

def run_command(command, description):
    """
    Executes a shell command and checks for errors.
    """
    print(f"--- Running {description} ---")
    # For commands that involve pipes, subprocess.run expects the command as a single string
    # and shell=True. This is generally less secure for arbitrary user input,
    # but for fixed internal commands like sed|awk, it's acceptable.
    # We'll print the command as a string for clarity.
    if isinstance(command, list):
        print(f"Command: {' '.join(command)}")
    else: # It's a string, likely for shell=True commands
        print(f"Command: {command}")

    try:
        # If the command is a string, assume shell=True is needed for pipes
        if isinstance(command, str):
            result = subprocess.run(command, check=True, capture_output=True, text=True, shell=True)
        else: # Otherwise, it's a list, run directly
            result = subprocess.run(command, check=True, capture_output=True, text=True)

        print("STDOUT:\n", result.stdout)
        if result.stderr:
            print("STDERR:\n", result.stderr)
        print(f"--- {description} completed successfully ---")
    except subprocess.CalledProcessError as e:
        print(f"Error during {description}:", file=sys.stderr)
        print(f"Command failed with exit code {e.returncode}", file=sys.stderr)
        print(f"STDOUT:\n{e.stdout}", file=sys.stderr)
        print(f"STDERR:\n{e.stderr}", file=sys.stderr)
        sys.exit(1)
    except FileNotFoundError:
        # This error is more likely for the first command in a pipe if shell=True is used
        # or for the direct command if shell=False.
        # We can't easily tell which command in a pipe failed if shell=True.
        print(f"Error: Command not found. Make sure necessary tools (sed, awk) are in your PATH.", file=sys.stderr)
        sys.exit(1)

def main():
    parser = argparse.ArgumentParser(
        description="Identify Reciprocal Best Hits (RBH) between two haploid genomes."
    )
    parser.add_argument(
        "--hapA_genome", required=True, help="Path to haploid genome A FASTA file."
    )
    parser.add_argument(
        "--hapB_genome", required=True, help="Path to haploid genome B FASTA file."
    )
    parser.add_argument(
        "--hapA_gff", required=True, help="Path to haploid genome A GFF file."
    )
    parser.add_argument(
        "--hapB_gff", required=True, help="Path to haploid genome B GFF file."
    )
    parser.add_argument(
        "-t", "--threads", type=int, default=1, help="Number of threads for BLAST."
    )
    parser.add_argument(
        "--evalue", type=float, default=1e-5, help="E-value cutoff for BLAST."
    )
    parser.add_argument(
        "--RBH_script", required=True, help="Path to the external RBH.py script."
    )
    parser.add_argument(
        "--out", required=True, help="Prefix for output RBH result file."
    )

    args = parser.parse_args()

    # Define intermediate filenames
    hapA_pep_fa = f"{os.path.splitext(args.hapA_genome)[0]}.pep.fa"
    hapB_pep_fa = f"{os.path.splitext(args.hapB_genome)[0]}.pep.fa"
    hapA_gene_fa = f"{os.path.splitext(args.hapA_genome)[0]}.gene.fa"
    hapB_gene_fa = f"{os.path.splitext(args.hapB_genome)[0]}.gene.fa"
    rbh_output_file = args.out # The user specified output name

    # --- Step 1: Extract protein sequences using agat ---
    print("\n##### Step 1: Extracting protein sequences with AGAT #####")

    # Command for haploid A
    agat_a_cmd = [
        "agat_sp_extract_sequences.pl",
        "-g", args.hapA_gff,
        "-f", args.hapA_genome,
        "-t", "cds",
        "-p",
        "-o", hapA_pep_fa
    ]
    run_command(agat_a_cmd, f"AGAT extraction for {args.hapA_genome}")

    # Command for haploid B
    agat_b_cmd = [
        "agat_sp_extract_sequences.pl",
        "-g", args.hapB_gff,
        "-f", args.hapB_genome,
        "-t", "cds",
        "-p",
        "-o", hapB_pep_fa
    ]
    run_command(agat_b_cmd, f"AGAT extraction for {args.hapB_genome}")

    # --- Step 1.5: Clean up FASTA headers using sed and awk ---
    print("\n##### Step 1.5: Cleaning up FASTA headers using sed and awk #####")

    # For haploid A
    # The command needs to be run with shell=True because of the pipes.
    # We construct the full shell command string.
    clean_a_cmd = f"sed 's/>.*gene=/>/g' {hapA_pep_fa} | sed 's/\\*//g' | awk '{{print $1}}' > {hapA_gene_fa}"
    run_command(clean_a_cmd, f"Cleaning headers for {hapA_pep_fa}")

    # For haploid B
    clean_b_cmd = f"sed 's/>.*gene=/>/g' {hapB_pep_fa} | sed 's/\\*//g' | awk '{{print $1}}' > {hapB_gene_fa}"
    run_command(clean_b_cmd, f"Cleaning headers for {hapB_pep_fa}")

    # --- Step 2: Call RBH.py ---
    print("\n##### Step 2: Calling RBH.py #####")

    rbh_cmd = [
        "python", args.RBH_script,
        "--hapA", hapA_gene_fa,
        "--hapB", hapB_gene_fa,
        "--out", rbh_output_file,
        "-t", str(args.threads),
        "--evalue", str(args.evalue)
    ]
    run_command(rbh_cmd, "RBH calculation")

    print(f"\n--- All steps completed successfully. RBH results saved to {rbh_output_file} ---")

if __name__ == "__main__":
    main()
