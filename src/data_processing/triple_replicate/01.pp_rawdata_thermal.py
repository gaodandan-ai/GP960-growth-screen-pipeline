# -*- coding: utf-8 -*-
"""
突变菌株胁迫条件，数据处理脚本，三个平行
处理突变菌株【基因过表达/TF失活/其他突变菌株】和正常对照菌株在两种条件conditon1和condition2【正常条件和胁迫条件下】的实验数据处理,生成四个CSV文件:

输入：
1. 原始实验数据文件夹，包含plate文件（Excel格式）
2. plate-gene映射文件（Excel格式）
3. 输出结果文件夹路径
4. condition1标签（nonstress）
5. condition2标签（stress）

输出：
1. mutant_{c1_label}_OD.csv - 突变菌株condition1条件下的OD数据
2. mutant_{c2_label}_OD.csv - 突变菌株condition2条件下的OD数据  
3. Ctrol_{c1_label}_OD.csv - 空质粒对照condition1条件下的OD数据
4. Ctrol_{c2_label}_OD.csv - 空质粒对照condition2条件下的OD数据



实验设计：
1 个 Excel 文件 = 1 个条件（Control = 正常，Treat = 高温）
每个 Excel 里包含 2 个原始板的数据：
1_* 开头的 sheet → 单数 Z 板（如 Z1）
2_* 开头的 sheet → 偶数 Z 板（如 Z2）
单数板和偶数板，都同时接种了两个条件：
单数板数据 → 存在于 Control.xlsx（正常）和 Treat.xlsx（高温）
偶数板数据 → 同样存在于 Control.xlsx 和 Treat.xlsx
H 行是单数板 + 偶数板的合并数据（一个 H sheet 包含两个原始板的 H 行）

具体来说：
- 每2个原始Z plate（单数与偶数）通过14个MTP实现，最后Condition2条件的文件为*_Treat.xlsx,Condition1条件的文件为*_Ctrl.xlsx，两个excel内sheet命名相同，意思为：
- 其中A-G行1-12列通过12个MTP实现（实验条件6，对照条件6），H行1-12列通过2个MTP实现（实验条件1，对照条件1）
-1_1_4: 原始单数Z plate A-G  1-4列，每列有3个重复，一共12列
-1_5_8: 原始单数Z plate A-G 5-8列，每列有3个重复，一共12列
-1_9_12: 原始单数Z plate A-G 9-12列，每列有3个重复，一共12列
-2_1_4: 原始偶数Z plate A-G 1-4列，每列有3个重复，一共12列
-2_5_8: 原始偶数Z plate A-G 5-8列，每列有3个重复，一共12列
-2_9_12: 原始偶数Z plate A-G 9-12列，每列有3个重复，一共12列

原始plate H行： 单数Z plate接种于右侧，偶数Z plate接种于左侧，一共有实验组和对照组两板，对应*_Treat.xlsx和*_Ctrl.xlsx的H sheet
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
        偶数Z plate 
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

对照组设计逻辑：（所有sheet相同）
Control板 *_Control.xlsx  Condition1
H1,H2,H3/H7,H8,H9为对照菌株（一组三个重复）
H4,H5,H6/H10,H11,H12为空白对照（一组三个重复）
Treat板 *_Treat.xlsx      Condition2
H1,H2,H3/H7,H8,H9为实验菌株（一组三个重复）
H4,H5,H6/H10,H11,H12为空白实验（一组三个重复）
"""


# -*- coding: utf-8 -*-
"""
突变菌株高温胁迫处理（匹配z1z2_Control/Treat.xlsx数据结构）
- 1个Excel=1个条件（Control=正常，Treat=高温）
- 每个Excel含7个sheet：1_*（单数板）、2_*（偶数板）、H（合并H行）
- 单数/偶数板均包含在两个条件的Excel中
"""

import pandas as pd
import numpy as np
from pathlib import Path

def load_plate_gene_mapping(mapping_file):
    """加载plate-gene映射表"""
    df = pd.read_excel(mapping_file)
    position_to_gene = {}
    for _, row in df.iterrows():
        # 映射表的location格式应为“板号+行+列”（如1A1、2H3）
        position = str(row['location']).strip()
        gene = str(row['gene']).strip()
        position_to_gene[position] = gene
    return position_to_gene, df

def extract_gene_data_from_sheet(df_sheet, position_to_gene, sheet_name):
    """
    从单个sheet提取基因数据（区分1_*/2_* sheet）
    sheet_name: 如1_1_4（单数板）、2_5_8（偶数板）、H（合并H行）
    """
    if df_sheet.shape[1] < 2:
        print(f"警告：{sheet_name} 列数不足，跳过")
        return {}, df_sheet.iloc[:, 1].dropna() if df_sheet.shape[1] > 0 else []

    # 提取时间列（第二列，第一列为序号）
    time_col = df_sheet.iloc[:, 1].dropna()
    max_rows = len(time_col)
    gene_data = {}

    # 1. 处理1_*/2_* sheet（A-G行）
    if sheet_name.startswith(('1_', '2_')):
        # 解析sheet名称：板类型（1=单数，2=偶数）、列范围（如1_4=1-4列）
        plate_type, start_col, end_col = sheet_name.split('_')
        plate_type = plate_type  # 1或2，用于拼接原始位置（如1A1、2B3）
        start_col = int(start_col)
        end_col = int(end_col)

        # 遍历A-G行
        for row_letter in ['A', 'B', 'C', 'D', 'E', 'F', 'G']:
            # 遍历列范围（如1-4列）
            for col in range(start_col, end_col + 1):
                # 原始板位置（如1A1、2C5）
                original_pos = f"{plate_type}{row_letter}{col}"
                if original_pos not in position_to_gene:
                    continue  # 映射表中没有的位置跳过
                gene = position_to_gene[original_pos]

                # 计算列偏移（如1-4列的偏移为0-3）
                col_offset = col - start_col
                # 每个基因3个重复（如偏移0→1-3列，偏移1→4-6列）
                for rep in range(3):
                    data_col = f"{row_letter}{col_offset * 3 + rep + 1}"
                    if data_col not in df_sheet.columns:
                        continue
                    # 提取数据并确保长度与时间列一致
                    data = df_sheet[data_col].dropna()
                    if len(data) == max_rows:
                        gene_data[f"{gene}-{rep+1}"] = data.values

    # 2. 处理H sheet（合并单数+偶数板的H行）
    elif sheet_name == 'H':
        # H行映射（单数板→1-6列，偶数板→7-12列）
        h_mapping = {
            # 单数板H行（对应A-F列1-6）
            '1': {
                1: ('A', [1,2,3]), 2: ('B', [1,2,3]), 3: ('C', [1,2,3]),
                4: ('D', [1,2,3]), 5: ('E', [1,2,3]), 6: ('F', [1,2,3]),
                7: ('A', [4,5,6]), 8: ('B', [4,5,6]), 9: ('C', [4,5,6]),
                10: ('D', [4,5,6]), 11: ('E', [4,5,6]), 12: ('F', [4,5,6])
            },
            # 偶数板H行（对应A-F列7-12）
            '2': {
                1: ('A', [7,8,9]), 2: ('B', [7,8,9]), 3: ('C', [7,8,9]),
                4: ('D', [7,8,9]), 5: ('E', [7,8,9]), 6: ('F', [7,8,9]),
                7: ('A', [10,11,12]), 8: ('B', [10,11,12]), 9: ('C', [10,11,12]),
                10: ('D', [10,11,12]), 11: ('E', [10,11,12]), 12: ('F', [10,11,12])
            }
        }

        # 遍历单数板（1）和偶数板（2）的H1-H12
        for plate_type in ['1', '2']:
            for h_pos in range(1, 13):  # H1到H12
                original_pos = f"{plate_type}H{h_pos}"
                if original_pos not in position_to_gene:
                    continue
                gene = position_to_gene[original_pos]

                # 获取对应的数据列（如1H1→A1/A2/A3）
                row_letter, cols = h_mapping[plate_type][h_pos]
                for rep, col in enumerate(cols):
                    data_col = f"{row_letter}{col}"
                    if data_col not in df_sheet.columns:
                        continue
                    data = df_sheet[data_col].dropna()
                    if len(data) == max_rows:
                        gene_data[f"{gene}-{rep+1}"] = data.values

    return gene_data, time_col

def extract_control_data_from_sheet(df_sheet, sheet_name):
    """
    从所有sheet提取对照数据
    - 对照菌株：H1/H2/H3/H7/H8/H9（所有sheet的H行）
    - 空白对照：H4/H5/H6/H10/H11/H12（所有sheet的H行）
    """
    if df_sheet.shape[1] < 2 or sheet_name != 'H':
        return {}, []  # 只从H sheet提取对照

    time_col = df_sheet.iloc[:, 1].dropna()
    max_rows = len(time_col)
    control_data = {}

    # 对照菌株位置（H1-H3、H7-H9）
    control_pos = ['H1', 'H2', 'H3', 'H7', 'H8', 'H9']
    for i, pos in enumerate(control_pos):
        if pos in df_sheet.columns:
            data = df_sheet[pos].dropna()
            if len(data) == max_rows:
                control_data[f"control_{i+1}"] = data.values

    # 空白对照位置（H4-H6、H10-H12）
    blank_pos = ['H4', 'H5', 'H6', 'H10', 'H11', 'H12']
    for i, pos in enumerate(blank_pos):
        if pos in df_sheet.columns:
            data = df_sheet[pos].dropna()
            if len(data) == max_rows:
                control_data[f"blank_{i+1}"] = data.values

    return control_data, time_col

def process_single_excel(excel_path, position_to_gene, condition_type):
    """处理单个Excel文件（1个文件=1个条件）"""
    print(f"\n=== 处理文件：{excel_path.name}（{condition_type}）===")
    try:
        excel_file = pd.ExcelFile(excel_path)
        sheet_names = excel_file.sheet_names
        print(f"包含sheet：{sheet_names}")
    except Exception as e:
        print(f"读取Excel失败：{e}")
        return {}, {}, None

    all_gene_data = {}
    all_control_data = {}
    time_points = None
    target_sheets = ['1_1_4', '1_5_8', '1_9_12', '2_1_4', '2_5_8', '2_9_12', 'H']

    for sheet in target_sheets:
        if sheet not in sheet_names:
            print(f"跳过不存在的sheet：{sheet}")
            continue
        try:
            # 读取当前sheet（第一行为表头）
            df_sheet = pd.read_excel(excel_path, sheet_name=sheet, header=0)
            print(f"  处理{sheet}：数据形状{df_sheet.shape}")

            # 提取基因数据
            gene_data, t_col = extract_gene_data_from_sheet(df_sheet, position_to_gene, sheet)
            # 提取对照数据（只从H sheet提）
            control_data, _ = extract_control_data_from_sheet(df_sheet, sheet)

            # 合并数据
            all_gene_data.update(gene_data)
            all_control_data.update(control_data)
            # 保存时间列（只取第一个有效sheet的时间）
            if time_points is None and len(t_col) > 0:
                time_points = t_col

            print(f"  {sheet}提取到{len(gene_data)}个基因，{len(control_data)}个对照数据")
        except Exception as e:
            print(f"  处理{sheet}出错：{e}")

    return all_gene_data, all_control_data, time_points

def main(rawdata_dir, plate_file, result_dir, c1_label='nonstress', c2_label='stress'):
    """
    主函数
    rawdata_dir: 原始数据目录（含Control/Treat.xlsx）
    plate_file: plate-gene映射表路径
    result_dir: 结果输出目录
    """
    # 初始化路径
    rawdata_dir = Path(rawdata_dir)
    result_dir = Path(result_dir)
    plate_file = Path(plate_file)
    result_dir.mkdir(parents=True, exist_ok=True)

    # 1. 加载plate-gene映射表
    if not plate_file.exists():
        print(f"错误：映射表{plate_file}不存在")
        return
    position_to_gene, _ = load_plate_gene_mapping(plate_file)
    print(f"加载映射表：共{len(position_to_gene)}个位置-基因对应关系")

    # 2. 找到Control和Treat文件
    control_files = list(rawdata_dir.glob("*Control*.xlsx"))
    treat_files = list(rawdata_dir.glob("*Treat*.xlsx"))
    if len(control_files) == 0:
        print("错误：未找到Control.xlsx文件")
        return
    if len(treat_files) == 0:
        print("错误：未找到Treat.xlsx文件")
        return
    control_file = control_files[0]
    treat_file = treat_files[0]

    # 3. 处理Control文件（Condition1=正常）
    gene_c1, ctrl_c1, times_c1 = process_single_excel(control_file, position_to_gene, "Condition1")
    # 处理Treat文件（Condition2=高温）
    gene_c2, ctrl_c2, times_c2 = process_single_excel(treat_file, position_to_gene, "Condition2")

    # 4. 统一时间列（用Control的时间，避免不一致）
    time_points = times_c1 if len(times_c1) > 0 else times_c2
    if len(time_points) == 0:
        print("错误：未提取到时间列")
        return

    # 5. 保存4个CSV结果
    def save_csv(data, filename):
        if len(data) == 0:
            print(f"跳过空文件：{filename}")
            return
        df = pd.DataFrame(data, index=time_points)
        df.index.name = "Time"  # 时间列命名
        output_path = result_dir / filename
        df.to_csv(output_path, encoding='utf-8-sig')
        print(f"\n保存结果：{output_path}")
        print(f"数据维度：时间点{len(time_points)} × 样本{len(data)}")

    # 输出4个文件
    save_csv(gene_c1, f"mutant_{c1_label}_OD.csv")    # 突变株-正常
    save_csv(gene_c2, f"mutant_{c2_label}_OD.csv")    # 突变株-高温
    save_csv(ctrl_c1, f"Control_{c1_label}_OD.csv")   # 对照株-正常
    save_csv(ctrl_c2, f"Control_{c2_label}_OD.csv")   # 对照株-高温

    print("\n=== 所有处理完成 ===")

