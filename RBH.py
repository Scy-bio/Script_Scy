import subprocess
import os
import argparse

def run_command(command, description):
    """
    Executes a shell command and prints its output.
    Raises an exception if the command fails.
    """
    print(f"--- Running: {description} ---")
    print(f"Command: {' '.join(command)}")
    try:
        result = subprocess.run(command, check=True, capture_output=True, text=True)
        print("STDOUT:\n", result.stdout)
        if result.stderr:
            print("STDERR:\n", result.stderr)
        print(f"--- {description} completed successfully ---")
    except subprocess.CalledProcessError as e:
        print(f"Error during {description}:")
        print(f"Command: {e.cmd}")
        print(f"Return Code: {e.returncode}")
        print(f"STDOUT:\n{e.stdout}")
        print(f"STDERR:\n{e.stderr}")
        raise
    except FileNotFoundError:
        print(f"Error: Command '{command[0]}' not found. "
              "Please ensure it's installed and in your system's PATH.")
        raise

def run_command_with_output_to_file(command, output_file, description):
    """
    Executes a shell command, pipes its stdout to a file, and prints its output.
    Raises an exception if the command fails.
    """
    print(f"--- Running: {description} (output to {output_file}) ---")
    print(f"Command: {' '.join(command)}")
    try:
        with open(output_file, 'w') as f:
            result = subprocess.run(command, check=True, stdout=f, stderr=subprocess.PIPE, text=True)
        if result.stderr:
            print("STDERR:\n", result.stderr)
        print(f"--- {description} completed successfully ---")
    except subprocess.CalledProcessError as e:
        print(f"Error during {description}:")
        print(f"Command: {e.cmd}")
        print(f"Return Code: {e.returncode}")
        print(f"STDERR:\n{e.stderr}")
        raise
    except FileNotFoundError:
        print(f"Error: Command '{command[0]}' not found. "
              "Please ensure it's installed and in your system's PATH.")
        raise

def find_rbhs(score_a2b_file, score_b2a_file, output_file):
    """
    Identifies Reciprocal Best Hits (RBH) from two BLAST score files.
    """
    print("\n##### Finding Reciprocal Best Hits (RBH) #####")
    blast_a2b = {}
    blast_b2a = {}

    # Read HAPA2HAPB.score
    with open(score_a2b_file, 'r') as f:
        for line in f:
            query, target, score = line.strip().split('\t')
            score = float(score)
            # Store the best hit for each query (highest bit score)
            if query not in blast_a2b or score > blast_a2b[query][1]:
                blast_a2b[query] = (target, score)

    # Read HAPB2HAPA.score
    with open(score_b2a_file, 'r') as f:
        for line in f:
            query, target, score = line.strip().split('\t')
            score = float(score)
            # Store the best hit for each query (highest bit score)
            if query not in blast_b2a or score > blast_b2a[query][1]:
                blast_b2a[query] = (target, score)

    # Find RBHs
    rbh_pairs = []
    for query_a, (target_b, score_a_b) in blast_a2b.items():
        # Check if target_b's best hit in the other direction is query_a
        if target_b in blast_b2a and blast_b2a[target_b][0] == query_a:
            rbh_pairs.append((query_a, target_b))

    # Write RBHs to output file
    with open(output_file, 'w') as f:
        for pair in rbh_pairs:
            f.write(f"{pair[0]}\t{pair[1]}\n")

    print(f"--- RBH identification completed. {len(rbh_pairs)} pairs found. Output saved to {output_file} ---")


def main():
    parser = argparse.ArgumentParser(description="Run DIAMOND BLAST and identify Reciprocal Best Hits (RBH).")
    parser.add_argument("--hapA", required=True, help="Path to HapA protein FASTA file.")
    parser.add_argument("--hapB", required=True, help="Path to HapB protein FASTA file.")
    parser.add_argument("--out", default="allele2blast.final.out", help="Output file for RBH pairs.")
    parser.add_argument("-t", "--threads", type=int, default=10, help="Number of threads for DIAMOND BLAST. Default: 10")
    parser.add_argument("--evalue", type=float, default=1e-10, help="E-value threshold for DIAMOND BLAST. Default: 1e-10")
    args = parser.parse_args()

    # Define file paths based on arguments
    hap_a_pep_fa = args.hapA
    hap_b_pep_fa = args.hapB
    hap_a_db = "HAPA"
    hap_b_db = "HAPB"
    outblast_a2b = "HAPA2HAPB.outblast"
    outblast_b2a = "HAPB2HAPA.outblast"
    score_a2b = "HAPA2HAPB.score"
    score_b2a = "HAPB2HAPA.score"
    rbh_output = args.out
    num_threads = str(args.threads) # Convert to string for subprocess
    e_value_threshold = str(args.evalue) # Convert to string for subprocess

    # --- Step 1: Create DIAMOND databases ---
    print("\n##### Creating DIAMOND Databases #####")
    run_command(["diamond", "makedb", "--in", hap_a_pep_fa, "--db", hap_a_db],
                f"Creating database for {hap_a_pep_fa}")
    run_command(["diamond", "makedb", "--in", hap_b_pep_fa, "--db", hap_b_db],
                f"Creating database for {hap_b_pep_fa}")

    # --- Step 2: Perform DIAMOND BLAST searches ---
    print("\n##### Performing DIAMOND BLAST Searches #####")
    # Dynamically build blast_params using the new arguments
    blast_params = [
        "--threads", num_threads,
        "--outfmt", "6",
        "--evalue", e_value_threshold,
        "--max-target-seqs", "20"
    ]

    run_command(["diamond", "blastp", "--db", f"{hap_b_db}.dmnd", "--query", hap_a_pep_fa, "--out", outblast_a2b] + blast_params,
                f"BLASTing {hap_a_pep_fa} against {hap_b_db}")
    run_command(["diamond", "blastp", "--db", f"{hap_a_db}.dmnd", "--query", hap_b_pep_fa, "--out", outblast_b2a] + blast_params,
                f"BLASTing {hap_b_pep_fa} against {hap_a_db}")

    # --- Step 3: Extract relevant columns for RBH analysis ---
    print("\n##### Extracting Scores for RBH Analysis #####")
    run_command_with_output_to_file(["cut", "-f", "1,2,12", outblast_a2b], score_a2b,
                                      f"Cutting columns from {outblast_a2b}")
    run_command_with_output_to_file(["cut", "-f", "1,2,12", outblast_b2a], score_b2a,
                                      f"Cutting columns from {outblast_b2a}")

    # --- Step 4: Find RBHs using custom Python function ---
    find_rbhs(score_a2b, score_b2a, rbh_output)

    print(f"\n--- All steps completed! RBH list saved to {rbh_output} ---")

if __name__ == "__main__":
    main()
