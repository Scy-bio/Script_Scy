import argparse
import os

def parse_annotation_file(filepath):
    """
    解析基因注释文件并将其存储在字典中。
    字典的 key: query ID, value: (seed_ortholog, Preferred_name)
    """
    annotations = {}
    if not os.path.exists(filepath):
        print(f"Error: Annotation file not found at '{filepath}'")
        return annotations

    with open(filepath, 'r') as f:
        # 读取并处理头部行，查找列索引
        header_line = f.readline().strip()
        if not header_line: # 检查是否是空文件
            print(f"Warning: Annotation file '{filepath}' is empty or has no header.")
            return annotations
        
        header = header_line.split('\t')
        
        try:
            query_idx = header.index('query')
            seed_idx = header.index('seed_ortholog')
            pref_name_idx = header.index('Preferred_name')
        except ValueError as e:
            print(f"Error: Missing expected column in header of '{filepath}': {e}. Expected 'query', 'seed_ortholog', 'Preferred_name'.")
            return annotations

        for line in f:
            parts = line.strip().split('\t')
            if len(parts) > max(query_idx, seed_idx, pref_name_idx):
                query_id = parts[query_idx]
                seed_ortholog = parts[seed_idx]
                preferred_name = parts[pref_name_idx]
                annotations[query_id] = (seed_ortholog, preferred_name)
            else:
                print(f"Warning: Skipping malformed line in '{filepath}': '{line.strip()}' (not enough columns)")
    return annotations

def main():
    parser = argparse.ArgumentParser(
        description="根据对应关系文件合并两个基因注释文件。",
        formatter_class=argparse.RawTextHelpFormatter # 保持帮助信息中的换行
    )
    parser.add_argument(
        "--hapA", 
        required=True, 
        help="CsiA基因注释文件的路径 (例如: 02.CsiA.anno2gene)"
    )
    parser.add_argument(
        "--hapB", 
        required=True, 
        help="CsiB基因注释文件的路径 (例如: 02.CsiB.anno2gene)"
    )
    parser.add_argument(
        "--col_rusult", 
        required=True, 
        help="CsiA和CsiB基因ID对应关系文件的路径 (例如: 02.col_rusult.txt)\n文件应为两列，第一列是CsiA ID，第二列是CsiB ID。"
    )
    parser.add_argument(
        "--out", 
        required=True, 
        help="输出合并后结果文件的路径 (例如: 03.result)"
    )

    args = parser.parse_args()

    print(f"开始解析 CsiA 注释文件: '{args.hapA}'...")
    csia_annotations = parse_annotation_file(args.hapA)
    if not csia_annotations:
        print("CsiA 注释文件解析失败或为空，程序退出。")
        return
    print(f"成功加载 {len(csia_annotations)} 条 CsiA 注释。")

    print(f"开始解析 CsiB 注释文件: '{args.hapB}'...")
    csib_annotations = parse_annotation_file(args.hapB)
    if not csib_annotations:
        print("CsiB 注释文件解析失败或为空，程序退出。")
        return
    print(f"成功加载 {len(csib_annotations)} 条 CsiB 注释。")

    print(f"正在处理对应关系文件: '{args.col_rusult}' 并写入结果到: '{args.out}'...")
    
    # 检查输出文件路径的目录是否存在，不存在则创建
    output_dir = os.path.dirname(args.out)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)
        print(f"创建输出目录: '{output_dir}'")

    try:
        with open(args.col_rusult, 'r') as col_f, open(args.out, 'w') as out_f:
            # 写入输出文件的头部
            out_f.write("CsiA_query\tCsiA_seed_ortholog\tCsiA_Preferred_name\t")
            out_f.write("CsiB_query\tCsiB_seed_ortholog\tCsiB_Preferred_name\n")

            processed_lines = 0
            for line in col_f:
                parts = line.strip().split('\t')
                if len(parts) == 2:
                    csia_id = parts[0]
                    csib_id = parts[1]

                    # 从字典中获取注释信息，如果找不到则用 '-' 填充
                    csia_info = csia_annotations.get(csia_id, ('-', '-'))
                    csib_info = csib_annotations.get(csib_id, ('-', '-'))

                    out_f.write(
                        f"{csia_id}\t{csia_info[0]}\t{csia_info[1]}\t"
                        f"{csib_id}\t{csib_info[0]}\t{csib_info[1]}\n"
                    )
                    processed_lines += 1
                else:
                    print(f"Warning: 对应关系文件 '{args.col_rusult}' 中发现格式不正确的行，已跳过: '{line.strip()}'")

        print(f"合并完成！共处理了 {processed_lines} 条对应关系。")
        print(f"结果已保存到: '{args.out}'")

    except FileNotFoundError:
        print(f"Error: 对应关系文件未找到或无法访问: '{args.col_rusult}'")
    except IOError as e:
        print(f"Error: 写入输出文件时发生I/O错误: '{args.out}' - {e}")
    except Exception as e:
        print(f"An unexpected error occurred during processing: {e}")

if __name__ == "__main__":
    main()
