import os
import sys
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def add_allele_ids_to_table(input_table_path: str, output_table_path: str) -> dict:
    """
    Reads an input table, adds a formatted 'AlleleXXXXX' ID to each line,
    and writes to an output table. It also returns a dictionary mapping
    original IDs (from the first two columns of the input table) to the new Allele IDs.

    Args:
        input_table_path (str): Path to the input annotation table.
        output_table_path (str): Path to the output table with Allele IDs.

    Returns:
        dict: A dictionary mapping original IDs (from col1 and col2) to new Allele IDs.
    """
    print(f"--- Step 1: Adding Allele IDs to '{input_table_path}' and generating '{output_table_path}' ---")
    id_to_allele_map = {}
    try:
        with open(input_table_path, 'r', encoding='utf-8') as infile, \
             open(output_table_path, 'w', encoding='utf-8') as outfile:
            for i, line in enumerate(infile, 1):
                original_line = line.strip()
                allele_id = f"Allele{i:05d}"
                new_line = f"{original_line}\t{allele_id}"
                outfile.write(new_line + '\n')

                # Assuming the first two columns of the original line are IDs to map
                parts = original_line.split('\t')
                if len(parts) >= 1:
                    id_to_allele_map[parts[0]] = allele_id
                if len(parts) >= 2:
                    id_to_allele_map[parts[1]] = allele_id
        print(f"Successfully generated '{output_table_path}'.")
    except FileNotFoundError:
        print(f"Error: Input table file '{input_table_path}' not found.", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Error processing '{input_table_path}': {e}", file=sys.stderr)
        sys.exit(1)
    return id_to_allele_map

def get_ids_from_allele_table(allele_table_path: str) -> tuple[set, set]:
    """
    Reads the 05.allele.table file and extracts IDs from the first and second columns.
    This replaces the previous 'homo' file logic.

    Args:
        allele_table_path (str): Path to the '05.allele.table' file.

    Returns:
        tuple[set, set]: Two sets of IDs (hapA_ids, hapB_ids).
    """
    print(f"\n--- Step 2: Reading IDs from '{allele_table_path}' for filtering ---")
    hapA_ids = set()
    hapB_ids = set()
    try:
        with open(allele_table_path, 'r', encoding='utf-8') as f:
            for line in f:
                parts = line.strip().split('\t') # Split by tab as it's a table
                if len(parts) >= 1:
                    hapA_ids.add(parts[0])
                if len(parts) >= 2:
                    hapB_ids.add(parts[1])
        print(f"Successfully extracted HapA and HapB IDs from '{allele_table_path}'.")
    except FileNotFoundError:
        print(f"Error: Allele table file '{allele_table_path}' not found. Make sure it's generated first.", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Error processing '{allele_table_path}': {e}", file=sys.stderr)
        sys.exit(1)
    return hapA_ids, hapB_ids

def filter_and_combine_fasta(
    hapA_fasta_path: str,
    hapB_fasta_path: str,
    hapA_ids_to_keep: set,
    hapB_ids_to_keep: set
) -> list[SeqRecord]:
    """
    Filters sequences from two FASTA files based on provided ID sets and combines them.

    Args:
        hapA_fasta_path (str): Path to HapA FASTA file.
        hapB_fasta_path (str): Path to HapB FASTA file.
        hapA_ids_to_keep (set): Set of IDs to keep from HapA FASTA.
        hapB_ids_to_keep (set): Set of IDs to keep from HapB FASTA.

    Returns:
        list[SeqRecord]: A list of combined and filtered SeqRecord objects.
    """
    print(f"\n--- Step 3: Filtering and combining FASTA sequences from '{hapA_fasta_path}' and '{hapB_fasta_path}' ---")
    combined_sequences = []
    try:
        # Filter HapA sequences
        for record in SeqIO.parse(hapA_fasta_path, "fasta"):
            if record.id in hapA_ids_to_keep:
                combined_sequences.append(record)
        print(f"Filtered {len(combined_sequences)} sequences from '{hapA_fasta_path}'.")

        # Filter HapB sequences
        initial_combined_len = len(combined_sequences)
        for record in SeqIO.parse(hapB_fasta_path, "fasta"):
            if record.id in hapB_ids_to_keep:
                combined_sequences.append(record)
        print(f"Filtered {len(combined_sequences) - initial_combined_len} sequences from '{hapB_fasta_path}'.")

        if not combined_sequences:
            print("Warning: No matching sequences found after filtering.", file=sys.stderr)

    except FileNotFoundError:
        print(f"Error: FASTA file '{hapA_fasta_path}' or '{hapB_fasta_path}' not found.", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Error processing FASTA files: {e}", file=sys.stderr)
        sys.exit(1)
    return combined_sequences

def rename_fasta_ids(
    sequences: list[SeqRecord],
    id_to_allele_map: dict,
    output_fasta_path: str
) -> list[SeqRecord]:
    """
    Renames the IDs of SeqRecord objects based on a provided mapping and writes
    the renamed sequences to a new FASTA file.

    Args:
        sequences (list[SeqRecord]): List of SeqRecord objects to rename.
        id_to_allele_map (dict): Dictionary mapping original IDs to new Allele IDs.
        output_fasta_path (str): Path to the output FASTA file with renamed IDs.

    Returns:
        list[SeqRecord]: A list of SeqRecord objects with renamed IDs.
    """
    print(f"\n--- Step 4: Renaming sequence IDs and generating '{output_fasta_path}' ---")
    renamed_sequences = []
    for record in sequences:
        original_id = record.id
        new_allele_id = id_to_allele_map.get(original_id, original_id) # Use original ID if not in map

        # Create a new SeqRecord with the updated ID
        new_record = SeqRecord(
            record.seq,
            id=new_allele_id,
            name=new_allele_id,
            description="" # Clear description or modify as needed
        )
        renamed_sequences.append(new_record)

    try:
        SeqIO.write(renamed_sequences, output_fasta_path, "fasta")
        print(f"Successfully generated '{output_fasta_path}' with renamed IDs.")
    except Exception as e:
        print(f"Error writing to '{output_fasta_path}': {e}", file=sys.stderr)
        sys.exit(1)
    return renamed_sequences


def main():
    """
    Main function to parse arguments and orchestrate the bioinformatics workflow.
    """
    parser = argparse.ArgumentParser(
        description="A script to process allele annotation tables and FASTA sequences."
    )
    parser.add_argument(
        '--result_file',
        type=str,
        required=True,
        help="Path to the input annotation table (e.g., 04.allele.anno.table)."
    )
    parser.add_argument(
        '--hapA_fasta',
        type=str,
        required=True,
        help="Path to the HapA gene FASTA file (e.g., ../03.allele.test/CskiameaA.gene.fa)."
    )
    parser.add_argument(
        '--hapB_fasta',
        type=str,
        required=True,
        help="Path to the HapB gene FASTA file (e.g., ../03.allele.test/CskiameaB.gene.fa)."
    )
    parser.add_argument(
        '--output_table',
        type=str,
        required=True,
        help="Path for the output table with Allele IDs (e.g., 05.allele.table)."
    )
    parser.add_argument(
        '--output_seq',
        type=str,
        required=True,
        help="Path for the output FASTA sequence file with renamed IDs (e.g., 05.allele.pep.fa.change)."
    )

    args = parser.parse_args()

    print("--- Starting prepare2needle.py Workflow ---")

    # Step 1: awk '{printf "%s\tAllele%05d\n", $0, NR}' 04.allele.anno.table > 05.allele.table
    # This also builds the ID mapping for later renaming
    id_to_allele_mapping = add_allele_ids_to_table(args.result_file, args.output_table)

    # Step 2: seqkit grep -f <(cut -f 1 05.allele.table) ... ; seqkit grep -f <(cut -f 2 05.allele.table) ...
    # Now, the IDs for grep come from the newly generated 05.allele.table
    hapA_ids_for_grep, hapB_ids_for_grep = get_ids_from_allele_table(args.output_table)

    # Step 3: cat CskiameaA.tmp1.fa CskiameaB.tmp1.fa > tmp.fa
    combined_fasta_records = filter_and_combine_fasta(
        args.hapA_fasta,
        args.hapB_fasta,
        hapA_ids_for_grep,
        hapB_ids_for_grep
    )

    # Step 4: awk '{print $1"\t"$3}' 05.allele.table > tmp1
    #         awk '{print $2"\t"$3}' 05.allele.table > tmp2
    #         cat tmp1 tmp2 > replace.id
    #         seqkit replace -p '(.+)' -r '{kv}' -k replace.id tmp.fa > 05.allele.pep.fa.change
    # This step uses the mapping from Step 1 and the combined records from Step 3
    rename_fasta_ids(
        combined_fasta_records,
        id_to_allele_mapping,
        args.output_seq
    )

    print("\n--- prepare2needle.py Workflow Completed Successfully ---")

if __name__ == "__main__":
    main()

