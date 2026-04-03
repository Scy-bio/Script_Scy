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
    print(f"Command: {' '.join(command)}")
    try:
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
        print(f"Error: Command not found. Make sure '{command[0]}' is in your PATH.", file=sys.stderr)
        sys.exit(1)

def main():
    parser = argparse.ArgumentParser(
        description="Automate WGDI analysis for synteny."
    )
    parser.add_argument(
        "--pA", required=True, help="Prefix for haploid A (e.g., hapA)."
    )
    parser.add_argument(
        "--pB", required=True, help="Prefix for haploid B (e.g., hapB)."
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
        "--conf_script", required=True, help="Path to the 01.generate_conf.py script."
    )
    parser.add_argument(
        "--conf_name", required=True, help="Base name for the WGDI config file (e.g., Csi will generate Csi.conf)."
    )
    parser.add_argument(
        "--blast_result", required=True, help="Path to the BLAST result file (e.g., HAPA2HAPB.outblast)."
    )
    parser.add_argument(
        "--genome1_name", required=True, help="Display name for genome 1 (e.g., C. sikamea(A))."
    )
    parser.add_argument(
        "--genome2_name", required=True, help="Display name for genome 2 (e.g., C. sikamea(B))."
    )
    parser.add_argument(
        "-t", "--threads", type=int, default=10, help="Number of threads/processes for WGDI collinearity (process parameter)."
    )
    parser.add_argument(
        "--evalue", type=float, default=1e-5, help="E-value cutoff for WGDI BLAST filter."
    )
    parser.add_argument(
        "--score", type=int, default=100, help="Score cutoff for WGDI BLAST filter."
    )


    args = parser.parse_args()

    # Define config file name
    conf_filename = f"{args.conf_name}.conf"

    # --- Step 1: Run 01.generate_conf.py to get necessary files ---
    print("\n##### Step 1: Running 01.generate_conf.py #####")

    # Command for haploid A
    gen_conf_a_cmd = [
        "python", args.conf_script,
        "-p", args.pA,
        args.hapA_genome,
        args.hapA_gff
    ]
    run_command(gen_conf_a_cmd, f"Generating WGDI files for {args.pA}")

    # Command for haploid B
    gen_conf_b_cmd = [
        "python", args.conf_script,
        "-p", args.pB,
        args.hapB_genome,
        args.hapB_gff
    ]
    run_command(gen_conf_b_cmd, f"Generating WGDI files for {args.pB}")

    # --- Step 2: Generate WGDI .conf file ---
    print("\n##### Step 2: Generating WGDI configuration file #####")

    config_content = f"""
[dotplot]
gff1 = {args.pA}.gff
gff2 = {args.pB}.gff
lens1 = {args.pA}.len
lens2 = {args.pB}.len
blast = {args.blast_result}
blast_reverse = false
genome1_name = {args.genome1_name}
genome2_name = {args.genome2_name}
multiple = 1
score = {args.score}
evalue = {args.evalue}
repeat_number = 10
position = order
ancestor_left = none
ancestor_top = none
markersize = 0.5
figsize = 10,10
savefig = wgdi.dot.pdf

[collinearity]
gff1 = {args.pA}.gff
gff2 = {args.pB}.gff
lens1 = {args.pA}.len
lens2 = {args.pB}.len
blast = {args.blast_result}
blast_reverse = false
multiple = 1
process = {args.threads}
evalue = {args.evalue}
score = {args.score}
grading = 50,40,25
mg = 40,40
pvalue = 0.2
repeat_number = 10
positon = order
savefile = wgdi.collinearity.txt
"""
    # Write the configuration to the specified file
    try:
        with open(conf_filename, 'w') as f:
            f.write(config_content.strip()) # .strip() removes leading/trailing whitespace
        print(f"WGDI configuration file '{conf_filename}' created successfully.")
    except IOError as e:
        print(f"Error writing WGDI config file '{conf_filename}': {e}", file=sys.stderr)
        sys.exit(1)

    # --- Step 3: Run WGDI ---
    print("\n##### Step 3: Running WGDI #####")

    # Run wgdi for dotplot
    wgdi_dot_cmd = ["wgdi", "-d", conf_filename]
    run_command(wgdi_dot_cmd, "WGDI dotplot analysis")

    # Run wgdi for collinearity
    wgdi_col_cmd = ["wgdi", "-icl", conf_filename]
    run_command(wgdi_col_cmd, "WGDI collinearity analysis")

    print(f"\n--- All WGDI steps completed successfully using '{conf_filename}' ---")

if __name__ == "__main__":
    main()
