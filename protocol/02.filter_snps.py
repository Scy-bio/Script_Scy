import pandas as pd
import argparse # Import the argparse module

def process_vcf_like_file(input_file, output_file, dp_threshold, gp_threshold, list_file=None):
    """
    Processes a VCF-like file, filters variants, generates statistics, and can output a site list.

    Args:
        input_file (str): Path to the input file.
        output_file (str): Path for the filtered output file.
        dp_threshold (int): Depth (DP) threshold.
        gp_threshold (float): Genotype call rate (GP) threshold (between 0 and 1).
        list_file (str, optional): Path for the output site list file (chr:pos format).
    """

    # Read the file. To handle an unknown number of columns, we'll read it line by line first.
    with open(input_file, 'r') as f:
        lines = f.readlines()

    if not lines:
        print("Input file is empty, cannot process.")
        return

    # Determine the number of individuals (tissues in this context)
    first_line_parts = lines[0].strip().split('\t')
    num_individuals = (len(first_line_parts) - 5) // 2
    print(f"Detected {num_individuals} tissues/samples in the file.")

    # Define column names
    column_names = ['Chromosome', 'Position', 'ID', 'REF', 'ALT']
    for i in range(num_individuals):
        column_names.append(f'Genotype_Individual_{i+1}')
        column_names.append(f'Depth_Individual_{i+1}')

    # Read data using pandas
    df = pd.read_csv(input_file, sep='\t', header=None, names=column_names)

    # Record the initial number of SNPs
    initial_snp_count = len(df)
    print(f"Initial SNP count: {initial_snp_count}")

    # Filtering logic
    filtered_rows = []
    site_list_entries = [] # To store "chr:pos" for the list file
    
    # 按照文章逻辑：至少有 20% 的组织样本深度达标
    min_dp_pass_rate = 0.20 

    for index, row in df.iterrows():
        called_gt_count = 0  # 记录成功检出基因型的组织数量
        dp_pass_count = 0    # 记录既成功检出，且深度达标(>= DP)的组织数量

        for i in range(num_individuals):
            genotype_col = f'Genotype_Individual_{i+1}'
            depth_col = f'Depth_Individual_{i+1}'

            # 只有当基因型非缺失时，才进行统计
            if row[genotype_col] != './.':
                called_gt_count += 1
                
                # 在检出基因型的前提下，检查深度是否达标
                try:
                    depth = int(row[depth_col])
                    if depth >= dp_threshold:
                        dp_pass_count += 1
                except ValueError:
                    pass

        # 计算整体检出率 (Overall CR)
        genotype_detection_rate = called_gt_count / num_individuals
        
        # 计算深度达标率 (有多少比例的组织其深度 >= dp_threshold)
        dp_pass_rate = dp_pass_count / num_individuals

        # 核心双重过滤逻辑：
        # 1. 整体检出率 >= gp_threshold (默认 0.5 即 50%)
        # 2. 深度达标率 >= min_dp_pass_rate (设置的 0.20 即 20%)
        if genotype_detection_rate >= gp_threshold and dp_pass_rate >= min_dp_pass_rate:
            filtered_rows.append(row)
            site_list_entries.append(f"{row['Chromosome']}:{row['Position']}")

    # Create the filtered DataFrame
    filtered_df = pd.DataFrame(filtered_rows, columns=column_names)

    # Record the number of SNPs after filtering
    filtered_snp_count = len(filtered_df)
    print(f"Filtered SNP count: {filtered_snp_count}")

    # Write the filtered data to a new file
    filtered_df.to_csv(output_file, sep='\t', header=False, index=False)
    print(f"Filtered data saved to: {output_file}")

    # Write the site list file if requested
    if list_file:
        with open(list_file, 'w') as f_list:
            for entry in site_list_entries:
                f_list.write(f"{entry}\n")
        print(f"Site list saved to: {list_file}")

    print("\n--- Statistics ---")
    print(f"Total number of tissues/samples: {num_individuals}")
    print(f"Total input SNPs: {initial_snp_count}")
    print(f"Total filtered SNPs: {filtered_snp_count}")
    print(f"Used sequencing depth threshold (DP): >={dp_threshold}")
    print(f"Used genotype call rate threshold (Overall CR): >={gp_threshold * 100:.1f}%")
    print(f"Used tissue depth pass rate threshold: >={min_dp_pass_rate * 100:.1f}%")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Filter SNPs in a VCF-like file based on strict multi-tissue thresholds and optionally output a site list.')

    parser.add_argument('--in', dest='input_file', required=True,
                        help='Path to the input file, e.g., T1.output.txt')

    parser.add_argument('--out', dest='output_file', required=True,
                        help='Path for the filtered output file, e.g., filtered_output.txt')

    parser.add_argument('--DP', dest='dp_threshold', type=int, default=10,
                        help=('Minimum sequencing depth per tissue (default is 10). '
                              'Note: The script strictly requires that at least 20%% of all sampled tissues '
                              'meet or exceed this depth threshold for a SNP to be retained.'))

    parser.add_argument('--GP', dest='gp_threshold', type=float, default=0.5,
                        help='Overall Genotype call rate threshold (default is 0.5, meaning call rate >= 50%%).')

    parser.add_argument('--list', dest='list_file', type=str, default=None,
                        help=('Optional: Path for the output site list file (chr:pos format). '
                              'e.g., filtered_sites.list. If not provided, no site list will be generated.'))

    # 解析命令行参数
    args = parser.parse_args()

    # 调用处理函数（只保留这一组调用即可）
    process_vcf_like_file(args.input_file, args.output_file, args.dp_threshold, args.gp_threshold, args.list_file)
