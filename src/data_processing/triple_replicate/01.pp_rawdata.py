# -*- coding: utf-8 -*-
"""
突变菌株胁迫条件，数据处理脚本，三个平行
处理突变菌株【基因过表达/TF失活/其他突变菌株】和正常对照菌株在两种条件conditon1和condition2【胁迫条件和正常条件下】的实验数据处理,生成四个CSV文件:

输入：
1. 原始实验数据文件夹，包含plate文件（Excel格式）
2. plate-gene映射文件（Excel格式）
3. 输出结果文件夹路径
4. condition1标签（例如: stress 或 nonstress）
5. condition2标签（例如: stress 或 nonstress）

输出：
1. mutant_{c1_label}_OD.csv - 突变菌株condition1条件下的OD数据
2. mutant_{c2_label}_OD.csv - 突变菌株condition2条件下的OD数据  
3. Ctrol_{c1_label}_OD.csv - 空质粒对照condition1条件下的OD数据
4. Ctrol_{c2_label}_OD.csv - 空质粒对照condition2条件下的OD数据



实验设计：
- 一共有34个plate.用于3000多个基因菌液的种子液培养。
- 每个原始plate通过7个MTP实现:MTP1(1-2), MTP2(3-4), MTP3(5-6), MTP4(7-8),MTP5(9-10),MTP6(11-12),MTP7(last)
- MTP1: 原始plate A-G 1-2列 -> A-G行1-6列(condition1),其中1-3列为原始plate中第1列的的三个重复,4-6列为原始plate中第2列的三个重复; 
                            A-G行7-12列(condition2),其中7-9列为原始plate中第1列的的三个重复,10-12列为原始plate中第2列的三个重复.
- MTP2: 原始plate A-G 3-4列 -> A-G行1-6列(condition1),其中1-3列为原始plate中第3列的的三个重复,4-6列为原始plate中第4列的三个重复; 
                            A-G行7-12列(condition2),其中7-9列为原始plate中第3列的的三个重复,10-12列为原始plate中第4列的三个重复.
- MTP3: 原始plate A-G 5-6列 -> A-G行1-6列(condition1),其中1-3列为原始plate中第5列的的三个重复,4-6列为原始plate中第6列的三个重复; 
                            A-G行7-12列(condition2),其中7-9列为原始plate中第5列的的三个重复,10-12列为原始plate中第6列的三个重复.
- MTP4: 原始plate A-G 7-8列 -> A-G行1-6列(condition1),其中1-3列为原始plate中第7列的的三个重复,4-6列为原始plate中第8列的三个重复; 
                            A-G行7-12列(condition2),其中7-9列为原始plate中第7列的的三个重复,10-12列为原始plate中第8列的三个重复.
- MTP5: 原始plate A-G 9-10列 -> A-G行1-6列(condition1),其中1-3列为原始plate中第9列的的三个重复,4-6列为原始plate中第10列的三个重复; 
                            A-G行7-12列(condition2),其中7-9列为原始plate中第9列的的三个重复,10-12列为原始plate中第10列的三个重复.
- MTP6: 原始plate A-G 11-12列 -> A-G行1-6列(condition1),其中1-3列为原始plate中第11列的的三个重复,4-6列为原始plate中第12列的三个重复; 
                            A-G行7-12列(condition2),其中7-9列为原始plate中第11列的的三个重复,10-12列为原始plate中第12列的三个重复.

- MTP7: 原始plate H行 1-12列 -> A-F行1-6列对应原始H1-12(condition1), 其中1-3列为原始plate中H行1-6列的三个重复,4-6列为原始plate中H行7-12列的三个重复;
                            A-F行7-12列对应原始H1-12(condition2),其中7-9列为原始plate中H行1-6列的三个重复,10-12列为原始plate中H行7-12列的三个重复.
        Noncondition2                           
        H1->A1,A2,A3
        H2->B1,B2,B3
        H3->C1,C2,C3
        H4->D1,D2,D3
        H5->E1,E2,E3
        H6->F1,F2,F3
        H7->A4,A5,A6
        H8->B4,B5,B6
        H9->C4,C5,C6
        H10->D4,D5,D6
        H11->E4,E5,E6
        H12->F4,F5,F6
        condition2 
        H1->A7,A8,A9
        H2->B7,B8,B9
        H3->C7,C8,C9
        H4->D7,D8,D9
        H5->E7,E8,E9
        H6->F7,F8,F9
        H7->A10,A11,A12
        H8->B10,B11,B12
        H9->C10,C11,C12
        H10->D10,D11,D12
        H11->E10,E11,E12
        H12->F10,F11,F12

对照组设计逻辑：
- 每个原始plate通过7个MTP实现:MTP1(1-2), MTP2(3-4), MTP3(5-6), MTP4(7-8),MTP5(9-10),MTP6(11-12),MTP7(last)
对于每个原始plate
原始plate A-G行1-2列 :condition1对照组是MTP1的H1-3,condition2对照组是MTP1的H7-9
原始plate A-G行3-4列 :condition1对照组是MTP2的H1-3,condition2对照组是MTP2的H7-9
原始plate A-G行5-6列 :condition1对照组是MTP3的H1-3,condition2对照组是MTP3的H7-9
原始plate A-G行7-8列 :condition1对照组是MTP4的H1-3,condition2对照组是MTP4的H7-9
原始plate A-G行9-10列 :condition1对照组是MTP5的H1-3,condition2对照组是MTP5的H7-9
原始plate A-G行11-12列 :condition1对照组是MTP6的H1-3,condition2对照组是MTP6的H7-9
原始plate H行1-12列 :condition1对照组是MTP7的H1-3,condition2对照组是MTP7的H7-9
"""

import pandas as pd
import numpy as np
import os
from pathlib import Path

def load_plate_gene_mapping(mapping_file):
    """加载plate-gene映射关系"""
    df = pd.read_excel(mapping_file)
    # 创建位置到基因的映射字典
    position_to_gene = {}
    for _, row in df.iterrows():
        position = row['location']
        gene = row['gene']
        position_to_gene[position] = gene
    return position_to_gene, df

def extract_gene_data_from_mtp(df_mtp, position_to_gene, plate_num, mtp_type):
    """从单个MTP sheet中提取基因数据，每个基因有三个重复"""
    # 检查DataFrame是否有足够的列
    if df_mtp.shape[1] < 2:
        print(f"警告: {mtp_type} sheet 列数不足，跳过处理")
        return {}, {}, df_mtp.iloc[:, 0].dropna() if df_mtp.shape[1] > 0 else []
    
    # 获取时间列
    time_col = df_mtp.iloc[:, 1].dropna()
    
    # 确保所有数据列都与时间列长度一致
    max_rows = len(time_col)
    
    gene_condition2_data = {}
    gene_condition1_data = {}
    
    if mtp_type in ['MTP1', 'MTP2', 'MTP3', 'MTP4', 'MTP5', 'MTP6']:
        # 计算原始plate中对应的列
        mtp_num = int(mtp_type[-1])  # 获取MTP编号
        col_offset = (mtp_num - 1) * 2  # MTP1:0, MTP2:2, MTP3:4, MTP4:6, MTP5:8, MTP6:10
        
        # 处理A-G行的基因数据
        for row_idx, row_letter in enumerate(['A', 'B', 'C', 'D', 'E', 'F', 'G']):
            for col_pair in range(2):  # 每个MTP处理2列
                # 原始plate中的位置
                original_col = col_pair + 1 + col_offset  # 计算原始列号
                original_pos = f"{int(plate_num)}{row_letter}{original_col}"
                
                if original_pos in position_to_gene:
                    gene = position_to_gene[original_pos]
                    
                    # condition1条件 - 每个基因有3个重复
                    for rep in range(3):
                        condition1_col = f"{row_letter}{col_pair * 3 + rep + 1}"  # 1-3列或4-6列
                        if condition1_col in df_mtp.columns:
                            condition1_data = df_mtp[condition1_col].dropna()
                            # 确保数据长度与时间列一致
                            if len(condition1_data) > max_rows:
                                condition1_data = condition1_data.iloc[:max_rows]
                            elif len(condition1_data) < max_rows:
                                # 如果数据长度不足，用最后一个值填充
                                last_val = condition1_data.iloc[-1] if len(condition1_data) > 0 else 0
                                padding = [last_val] * (max_rows - len(condition1_data))
                                condition1_data = pd.concat([condition1_data, pd.Series(padding)])
                            gene_rep_name = f"{gene}-{rep + 1}"
                            gene_condition1_data[gene_rep_name] = condition1_data.values
                    
                    # condition2条件 - 每个基因有3个重复
                    for rep in range(3):
                        condition2_col = f"{row_letter}{col_pair * 3 + rep + 7}"  # 7-9列或10-12列
                        if condition2_col in df_mtp.columns:
                            condition2_data = df_mtp[condition2_col].dropna()
                            # 确保数据长度与时间列一致
                            if len(condition2_data) > max_rows:
                                condition2_data = condition2_data.iloc[:max_rows]
                            elif len(condition2_data) < max_rows:
                                # 如果数据长度不足，用最后一个值填充
                                last_val = condition2_data.iloc[-1] if len(condition2_data) > 0 else 0
                                padding = [last_val] * (max_rows - len(condition2_data))
                                condition2_data = pd.concat([condition2_data, pd.Series(padding)])
                            gene_rep_name = f"{gene}-{rep + 1}"
                            gene_condition2_data[gene_rep_name] = condition2_data.values
    
    elif mtp_type == 'MTP7':
        # MTP7处理原始plate的H行
        # H1-H12的映射关系
        h_mapping = {
            # condition1 mapping
            1: ('A', [1, 2, 3]), 2: ('B', [1, 2, 3]), 3: ('C', [1, 2, 3]), 
            4: ('D', [1, 2, 3]), 5: ('E', [1, 2, 3]), 6: ('F', [1, 2, 3]),
            7: ('A', [4, 5, 6]), 8: ('B', [4, 5, 6]), 9: ('C', [4, 5, 6]),
            10: ('D', [4, 5, 6]), 11: ('E', [4, 5, 6]), 12: ('F', [4, 5, 6])
        }
        
        for h_idx in range(12):  # H1-H12
            original_pos = f"{int(plate_num)}H{h_idx + 1}"
            if original_pos in position_to_gene:
                gene = position_to_gene[original_pos]
                
                # 获取映射信息
                row_letter, condition1_cols = h_mapping[h_idx + 1]
                
                # condition1条件 - 三个重复
                for rep in range(3):
                    condition1_col = f"{row_letter}{condition1_cols[rep]}"
                    if condition1_col in df_mtp.columns:
                        condition1_data = df_mtp[condition1_col].dropna()
                        # 确保数据长度与时间列一致
                        if len(condition1_data) > max_rows:
                            condition1_data = condition1_data.iloc[:max_rows]
                        elif len(condition1_data) < max_rows:
                            # 如果数据长度不足，用最后一个值填充
                            last_val = condition1_data.iloc[-1] if len(condition1_data) > 0 else 0
                            padding = [last_val] * (max_rows - len(condition1_data))
                            condition1_data = pd.concat([condition1_data, pd.Series(padding)])
                        gene_rep_name = f"{gene}-{rep + 1}"
                        gene_condition1_data[gene_rep_name] = condition1_data.values
                
                # condition2条件 - 三个重复
                condition2_cols = [col + 6 for col in condition1_cols]  # 7-9列或10-12列
                for rep in range(3):
                    condition2_col = f"{row_letter}{condition2_cols[rep]}"
                    if condition2_col in df_mtp.columns:
                        condition2_data = df_mtp[condition2_col].dropna()
                        # 确保数据长度与时间列一致
                        if len(condition2_data) > max_rows:
                            condition2_data = condition2_data.iloc[:max_rows]
                        elif len(condition2_data) < max_rows:
                            # 如果数据长度不足，用最后一个值填充
                            last_val = condition2_data.iloc[-1] if len(condition2_data) > 0 else 0
                            padding = [last_val] * (max_rows - len(condition2_data))
                            condition2_data = pd.concat([condition2_data, pd.Series(padding)])
                        gene_rep_name = f"{gene}-{rep + 1}"
                        gene_condition2_data[gene_rep_name] = condition2_data.values
    
    return gene_condition2_data, gene_condition1_data, time_col

def extract_control_data_from_mtp(df_mtp, mtp_type, gene_mapping_df=None, plate_num=None):
    """从MTP中提取对照组数据
    
    对照组设计逻辑：
    - MTP1的对照逻辑：为A-G行1-2列的基因提供对照
      H1-3对应A-G行1-2列的基因（condition1条件）
      H7-9对应A-G行1-2列的基因（condition2条件）
    - MTP2的对照逻辑：为A-G行3-4列的基因提供对照
      H1-3对应A-G行3-4列的基因（condition1条件）
      H7-9对应A-G行3-4列的基因（condition2条件）
    - 依此类推到MTP6
    - MTP7的对照逻辑：为H行1-12列的基因提供对照
      H1-H3对应H行1-12列的基因（condition1条件）
      H7-H9对应H行1-12列的基因（condition2条件）
    """
    # 检查DataFrame是否有足够的列
    if df_mtp.shape[1] < 2:
        print(f"警告: {mtp_type} sheet 列数不足，跳过对照组处理")
        return {}, {}, df_mtp.iloc[:, 0].dropna() if df_mtp.shape[1] > 0 else []
    
    # 获取时间列
    time_col = df_mtp.iloc[:, 1].dropna()
    # 确保所有数据列都与时间列长度一致
    max_rows = len(time_col)
    
    control_condition2_data = {}
    control_condition1_data = {}
    
    # 获取该MTP对应的目标基因列表
    def get_target_genes_for_mtp(mtp_type, gene_mapping_df, plate_num):
        """获取该MTP对应的目标基因列表"""
        if gene_mapping_df is None:
            return []
        
        # 解析位置信息
        import re
        def parse_location(location):
            match = re.match(r'(\d+)([A-H])(\d+)', location)
            if match:
                return int(match.group(1)), match.group(2), int(match.group(3))
            return None, None, None
        
        # 创建副本避免修改原始数据
        gene_mapping_copy = gene_mapping_df.copy()
        gene_mapping_copy[['plate', 'row', 'col']] = gene_mapping_copy['location'].apply(parse_location).apply(pd.Series)
        
        # 将plate_num转换为整数进行比较
        plate_id = int(plate_num)
        # 只获取当前plate的基因
        current_plate_data = gene_mapping_copy[gene_mapping_copy['plate'] == plate_id]
        
        mtp_num = int(mtp_type.replace('MTP', ''))
        target_genes = []
        
        if mtp_num <= 6:  # MTP1-MTP6
            # 对于MTP1-6，收集该MTP对应列的所有基因
            target_cols = [2*mtp_num-1, 2*mtp_num]  # MTP1->[1,2], MTP2->[3,4], etc.
            
            for row in ['A', 'B', 'C', 'D', 'E', 'F', 'G']:
                for col in target_cols:
                    gene_data = current_plate_data[(current_plate_data['row'] == row) & (current_plate_data['col'] == col)]
                    if not gene_data.empty:
                        target_genes.extend(gene_data['gene'].tolist())
                        
        else:  # MTP7
            # 对于MTP7，收集H行1-12列的基因
            for col_idx in range(12):  # 1-12列
                gene_data = current_plate_data[(current_plate_data['row'] == 'H') & (current_plate_data['col'] == col_idx + 1)]
                if not gene_data.empty:
                    target_genes.extend(gene_data['gene'].tolist())
        
        return target_genes
    
    # 获取目标基因列表
    target_genes = get_target_genes_for_mtp(mtp_type, gene_mapping_df, plate_num)
    # 生成基因列表字符串（用分号分隔）
    gene_list = ';'.join(target_genes) if target_genes else "NoGenes"
    
    # 提取对照组数据：H1-H3为condition1，H7-H9为condition2
    for rep in range(3):
        # condition1对照组 (H1-H3)
        condition1_col = f"H{rep + 1}"
        if condition1_col in df_mtp.columns:
            condition1_data = df_mtp[condition1_col].dropna()
            # 确保数据长度与时间列一致
            if len(condition1_data) > max_rows:
                condition1_data = condition1_data.iloc[:max_rows]
            elif len(condition1_data) < max_rows:
                # 如果数据长度不足，用最后一个值填充
                last_val = condition1_data.iloc[-1] if len(condition1_data) > 0 else 0
                padding = [last_val] * (max_rows - len(condition1_data))
                condition1_data = pd.concat([condition1_data, pd.Series(padding)])
            #
            col_name = f"ctrol{rep + 1}-{gene_list}"
            control_condition1_data[col_name] = condition1_data.values
            print(f"    提取 condition1 对照: {condition1_col} -> {col_name[:50]}...")
        
        # condition2对照组 (H7-H9)
        condition2_col = f"H{rep + 7}"
        if condition2_col in df_mtp.columns:
            condition2_data = df_mtp[condition2_col].dropna()
            # 确保数据长度与时间列一致
            if len(condition2_data) > max_rows:
                condition2_data = condition2_data.iloc[:max_rows]
            elif len(condition2_data) < max_rows:
                # 如果数据长度不足，用最后一个值填充
                last_val = condition2_data.iloc[-1] if len(condition2_data) > 0 else 0
                padding = [last_val] * (max_rows - len(condition2_data))
                condition2_data = pd.concat([condition2_data, pd.Series(padding)])
            # 使用与pp_MSG.py一致的命名格式
            col_name = f"ctrol{rep + 1}-{gene_list}"
            control_condition2_data[col_name] = condition2_data.values
            print(f"    提取 condition2 对照: {condition2_col} -> {col_name[:50]}...")
    
    return control_condition2_data, control_condition1_data, time_col

def process_single_plate(file_path, position_to_gene, plate_num, gene_mapping_df=None):
    """处理单个plate文件"""
    print(f"处理文件: {file_path}")
    
    # 读取Excel文件的所有sheet
    try:
        excel_file = pd.ExcelFile(file_path)
        sheet_names = excel_file.sheet_names
        print(f"发现sheets: {sheet_names}")
    except Exception as e:
        print(f"读取文件失败: {e}")
        return {}, {}, {}, {}, None
    
    all_gene_condition2 = {}
    all_gene_condition1 = {}
    all_control_condition2 = {}
    all_control_condition1 = {}
    time_column = None
    
    # 处理7个MTP
    mtp_ranges = {
        1: '(1-2)',
        2: '(3-4)', 
        3: '(5-6)',
        4: '(7-8)',
        5: '(9-10)',
        6: '(11-12)',
        7: '(last)'
    }
    
    for mtp_num in range(1, 8):
        mtp_name = f'MTP{mtp_num}'
        target_range = mtp_ranges[mtp_num]
        
        # 寻找包含对应范围的sheet名
        actual_sheet_name = None
        for s in sheet_names:
            if target_range in s:
                actual_sheet_name = s
                break
        
        if actual_sheet_name:
            try:
                df_mtp = pd.read_excel(file_path, sheet_name=actual_sheet_name, header=0)
                print(f"处理 {actual_sheet_name} (作为 {mtp_name}), 形状: {df_mtp.shape}")
                
                # 提取基因数据
                gene_condition2, gene_condition1, time_col = extract_gene_data_from_mtp(
                    df_mtp, position_to_gene, plate_num, mtp_name
                )
                
                # 提取对照组数据，传递基因映射数据和plate编号
                control_condition2, control_condition1, _ = extract_control_data_from_mtp(
                    df_mtp, mtp_name, gene_mapping_df, plate_num
                )
                
                # 合并数据
                all_gene_condition2.update(gene_condition2)
                all_gene_condition1.update(gene_condition1)
                all_control_condition2.update(control_condition2)
                all_control_condition1.update(control_condition1)
                
                # 保存时间列
                if time_column is None and len(time_col) > 0:
                    time_column = time_col
                    
            except Exception as e:
                print(f"处理 {actual_sheet_name} 时出错: {e}")
        else:
            print(f"未找到sheet: {actual_sheet_name}")
    
    return all_gene_condition2, all_gene_condition1, all_control_condition2, all_control_condition1, time_column

def main():
    import argparse
    parser = argparse.ArgumentParser(description='处理高通量生长曲线数据')
    parser.add_argument('--rawdata_dir', type=str, required=True,
                        help='原始实验数据目录，如: /data/zuoll/1.project/02.GOE/GP960-growth-screen-pipeline/data/raw/single_replicate/highmethanol')
    parser.add_argument('--plate_file', type=str, required=True,
                        help='实验对应关系文件，如: /data/zuoll/1.project/02.GOE/GP960-growth-screen-pipeline/data/raw/single_replicate/plate-gene-mapping.xlsx')
    parser.add_argument('--result_dir', type=str, required=True,
                        help='处理后的结果目录，如: /data/zuoll/1.project/02.GOE/GP960-growth-screen-pipeline/data/results/single_replicate/highmethanol')
    parser.add_argument('--c1_label', type=str, default='stress',
                    choices=['stress', 'nonstress'],  # 限制只能选择这两个值
                    help='Condition 1 的标签: stress (胁迫) 或 nonstress (无胁迫)')

    parser.add_argument('--c2_label', type=str, default='nonstress',
                    choices=['stress', 'nonstress'],  # 同样限制选择范围
                    help='Condition 2 的标签: stress (胁迫) 或 nonstress (无胁迫)')
    args = parser.parse_args()
    
    # 设置路径
    rawdata_dir = Path(args.rawdata_dir)
    plate_file_path = Path(args.plate_file)
    result_dir = Path(args.result_dir)
    
    # 创建结果目录
    result_dir.mkdir(parents=True, exist_ok=True)
    
    # 加载plate-gene映射
    if not plate_file_path.exists():
        print(f"错误: 找不到 mapping 文件 {plate_file_path}")
        return

    position_to_gene, gene_mapping_df = load_plate_gene_mapping(plate_file_path)
    print(f"加载了 {len(position_to_gene)} 个基因映射关系")
    
    # 初始化数据收集器
    all_gene_condition2_data = {}
    all_gene_condition1_data = {}
    all_control_condition2_data = {}
    all_control_condition1_data = {}
    time_points = None
    
    # 自动识别plate文件
    # 查找符合模式 {plate_num:02d}-{args.condition}.xlsx 的文件
    # 或者如果不指定 condition 参数，查找所有 .xlsx 文件并尝试解析
    
    # 移除 condition 参数依赖，改用自动搜索
    xlsx_files = sorted(list(rawdata_dir.glob('*.xlsx')))
    
    # 过滤掉 mapping 文件
    xlsx_files = [f for f in xlsx_files if "mapping" not in f.name.lower()]
    
    if not xlsx_files:
        print(f"警告: 在 {rawdata_dir} 未找到数据文件")
        return

    print(f"找到 {len(xlsx_files)} 个数据文件: {[f.name for f in xlsx_files]}")
    
    # 处理找到的每个文件
    for plate_file in xlsx_files:
        # 尝试从文件名提取 plate_num (假设文件名以数字开头，如 01-xxx.xlsx)
        try:
            plate_prefix = plate_file.name.split('-')[0]
            plate_num = int(plate_prefix)
        except ValueError:
            print(f"警告: 无法从文件名 {plate_file.name} 提取 plate 编号，跳过处理")
            continue
            
        gene_condition2, gene_condition1, control_condition2, control_condition1, times = process_single_plate(
            plate_file, position_to_gene, plate_num, gene_mapping_df)
        
        # 合并数据
        all_gene_condition2_data.update(gene_condition2)
        all_gene_condition1_data.update(gene_condition1)
        all_control_condition2_data.update(control_condition2)
        all_control_condition1_data.update(control_condition1)
        
        if time_points is None:
            time_points = times
            
        print(f"从 {plate_file.name} 提取了 {len(gene_condition2)} 个condition2基因, {len(gene_condition1)} 个condition1 基因")
    
    # 创建DataFrame并保存
    print("生成CSV文件...")
    print(f"总共收集到: condition1基因 {len(all_gene_condition1_data)}, condition2 基因 {len(all_gene_condition2_data)}")
    print(f"对照数据: condition1 {len(all_control_condition1_data)}, condition2  {len(all_control_condition2_data)}")
    
    # 1. mutant_condition1_OD.csv (使用c1_label)
    if all_gene_condition1_data:
        df_gene_condition1 = pd.DataFrame(all_gene_condition1_data, index=time_points)
        filename = f"mutant_{args.c1_label}_OD.csv"
        df_gene_condition1.to_csv(result_dir / filename)
        print(f"生成 {filename}: {df_gene_condition1.shape}")
    
    # 2. mutant_condition2_OD.csv (使用c2_label)
    if all_gene_condition2_data:
        df_gene_condition2 = pd.DataFrame(all_gene_condition2_data, index=time_points)
        filename = f"mutant_{args.c2_label}_OD.csv"
        df_gene_condition2.to_csv(result_dir / filename) 
        print(f"生成 {filename}: {df_gene_condition2.shape}")
    
    # 3. Ctrol_condition1_OD.csv (使用c1_label)
    if all_control_condition1_data:
        df_control_condition1 = pd.DataFrame(all_control_condition1_data, index=time_points)
        filename = f"Ctrol_{args.c1_label}_OD.csv"
        df_control_condition1.to_csv(result_dir / filename)  
        print(f"生成 {filename}: {df_control_condition1.shape}")

    # 4. Ctrol_condition2_OD.csv (使用c2_label)
    if all_control_condition2_data:
        df_control_condition2 = pd.DataFrame(all_control_condition2_data, index=time_points)
        filename = f"Ctrol_{args.c2_label}_OD.csv"
        df_control_condition2.to_csv(result_dir / filename)
        print(f"生成 {filename}: {df_control_condition2.shape}")

if __name__ == "__main__":
    main()
