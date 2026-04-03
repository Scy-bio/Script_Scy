#!/usr/bin/env python3

import argparse
import pandas as pd
import re
import sys

def process_emapper_annotations(input_file_path, output_file_path):
    """
    Processes an eggNOG-mapper annotation file to extract specific columns:
    'query', 'seed_ortholog', and 'Preferred_name'.
    It handles inconsistent spacing and potential non-standard delimiters.

    Args:
        input_file_path (str): The path to the input .emapper.annotations file.
        output_file_path (str): The path where the extracted data will be saved.
    """
    try:
        data_rows = []
        header_found = False
        header = []
        
        # Define the columns we want to extract
        columns_to_extract = ['query', 'seed_ortholog', 'Preferred_name']

        with open(input_file_path, 'r', encoding='utf-8') as f_in:
            for line_num, line in enumerate(f_in, 1):
                # Skip comment lines (starting with '##')
                if line.startswith('##'):
                    continue

                # Identify the header line (starts with '#query')
                if line.startswith('#query'):
                    # Remove the '#' prefix, strip leading/trailing whitespace
                    # Replace common problematic spaces with standard spaces
                    clean_header_line = line[1:].strip().replace('\xa0', ' ').replace(' ', ' ')
                    
                    # Split header by a single tab, then clean each part by collapsing multiple spaces
                    header = [' '.join(h.split()).strip() for h in clean_header_line.split('\t')]
                    
                    # Remove any empty strings that might result from extra tabs in header
                    header = [h for h in header if h]
                    
                    if not header:
                        sys.stderr.write(f"Error: Header line found but no valid column names parsed from line {line_num}: '{line.strip()}'\n")
                        sys.exit(1)

                    header_found = True
                    continue

                # Process data lines after the header has been found
                if header_found:
                    # Clean the line by replacing problematic spaces
                    cleaned_line = line.strip().replace('\xa0', ' ').replace(' ', ' ')
                    
                    # If the line is empty after cleaning, skip it
                    if not cleaned_line:
                        continue

                    # Split by tab. We assume tabs are the primary column delimiters.
                    parts = cleaned_line.split('\t')

                    # Clean up each part (field) by collapsing multiple internal spaces and stripping
                    parts = [' '.join(p.split()).strip() for p in parts]
                    
                    # Remove any empty strings that might result from extra tabs or failed splits
                    parts = [p for p in parts if p]

                    # Validate the number of columns
                    if len(parts) > len(header):
                        # If more parts than header, it might be due to a tab inside a field.
                        # We truncate to match the expected number of columns.
                        sys.stderr.write(f"Warning: Line {line_num} has more columns ({len(parts)}) than header ({len(header)}). Truncating.\n")
                        # sys.stderr.write(f"Original Line: '{line.strip()}'\n") # Uncomment for detailed debugging
                        # sys.stderr.write(f"Parsed Parts: {parts}\n") # Uncomment for detailed debugging
                        parts = parts[:len(header)]
                    elif len(parts) < len(header):
                        # If fewer parts, it means some columns are missing or combined.
                        sys.stderr.write(f"Warning: Line {line_num} has fewer columns ({len(parts)}) than header ({len(header)}). Skipping line: '{line.strip()}'\n")
                        continue # Skip lines that don't match the expected column count

                    data_rows.append(parts)

        if not header_found:
            raise ValueError(f"Error: Could not find the header line starting with '#query' in '{input_file_path}'. "
                             "Please ensure it's a valid emapper annotation file.")
        
        if not data_rows:
            print(f"Warning: No valid data rows found after header in '{input_file_path}'. Output file will be empty.")
            # Create an empty DataFrame with expected columns to avoid errors if no data
            df = pd.DataFrame(columns=header)
        else:
            df = pd.DataFrame(data_rows, columns=header)

        # Final check if all required columns for extraction exist in the DataFrame
        missing_columns = [col for col in columns_to_extract if col not in df.columns]
        if missing_columns:
            sys.stderr.write(f"Error: The following required columns were not found in the parsed data: {', '.join(missing_columns)}\n")
            sys.stderr.write(f"Available columns in parsed data: {df.columns.tolist()}\n")
            raise ValueError("Crucial columns missing after parsing. Check input file format or column names.")

        # Select only the desired columns
        extracted_df = df[columns_to_extract]

        # Save the extracted data to the output file as tab-separated
        extracted_df.to_csv(output_file_path, sep='\t', index=False)
        print(f"Successfully extracted columns and saved to '{output_file_path}'")

    except FileNotFoundError:
        sys.stderr.write(f"Error: Input file '{input_file_path}' not found.\n")
        sys.exit(1)
    except Exception as e:
        sys.stderr.write(f"An unexpected error occurred: {e}\n")
        sys.exit(1)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Extracts 'query', 'seed_ortholog', and 'Preferred_name' columns "
                    "from an eggNOG-mapper annotation file, handling varied internal spacing."
    )
    parser.add_argument(
        "--in",
        dest="input_file",
        required=True,
        help="Path to the input eggNOG-mapper annotations file (e.g., 01.CsiA.emapper.annotations)."
    )
    parser.add_argument(
        "--out",
        dest="output_file",
        required=True,
        help="Path to the output file where the extracted data will be saved (e.g., 02.CsiA.anno2gene)."
    )

    args = parser.parse_args()

    process_emapper_annotations(args.input_file, args.output_file)
