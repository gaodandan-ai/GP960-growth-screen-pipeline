#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Thermal Overpress All 数据处理脚本
处理 Thermal_40 和 Control_30 实验数据, 生成四个 CSV 文件

输入：
1. ALL.xlsx (包含 4 个数据 sheet 和 1 个 meta sheet)
2. plate_gene_mapping.xlsx (映射数据)

输出：
1. 01.mutant_{c1_label}_OD.csv
2. 02.mutant_{c2_label}_OD.csv
3. 03.Ctrol_{c1_label}_OD.csv
4. 04.Ctrol_{c2_label}_OD.csv
"""

import pandas as pd
import numpy as np
import os
from pathlib import Path
import argparse

def load_plate_gene_mapping(mapping_file):
    """加载 plate-gene 映射关系"""
    print(f"加载映射文件: {mapping_file}")
    try:
        # 尝试读取 'meta' sheet，如果不存在则读取第一个 sheet
        xls = pd.ExcelFile(mapping_file)
        sheet_name = 'meta' if 'meta' in xls.sheet_names else xls.sheet_names[0]
        df = pd.read_excel(xls, sheet_name=sheet_name)
    except Exception as e:
        print(f"读取映射文件失败: {e}")
        return {}

    # 创建 Position 到 old_locus_tag 的映射
    # 注意：映射文件中列名为 'Position'，原始数据中列名也对应
    position_to_gene = {}
    if 'Position' in df.columns and 'old_locus_tag' in df.columns:
        for _, row in df.iterrows():
            pos = str(row['Position']).strip()
            gene = str(row['old_locus_tag']).strip()
            if pos and gene and pos != 'nan' and gene != 'nan':
                position_to_gene[pos] = gene
    else:
        print(f"警告: 映射文件缺失 'Position' 或 'old_locus_tag' 列. 现有列: {df.columns.tolist()}")
        
    return position_to_gene

def process_mutant_sheet(df, position_to_gene):
    """处理突变株数据 sheet"""
    if 'Time(hours)' not in df.columns:
        print("错误: 数据中缺少 'Time(hours)' 列")
        return None, None
    
    time_points = df['Time(hours)'].values
    mutant_data = {}
    
    # 遍历所有列，如果是位置则映射到基因名
    for col in df.columns:
        if col == 'Time(hours)':
            continue
        
        if col in position_to_gene:
            gene_name = position_to_gene[col]
            mutant_data[gene_name] = df[col].values
        else:
            # print(f"跳过未识别的位置: {col}")
            pass
            
    return mutant_data, time_points

def process_vc_sheet(df):
    """处理 VC 菌株数据 sheet，取均值"""
    if 'Time(hours)' not in df.columns:
        return None
    
    # 提取所有数据列（排除时间列和全是 NaN 的列）
    data_cols = [col for col in df.columns if col != 'Time(hours)']
    vc_df = df[data_cols].dropna(axis=1, how='all')
    
    # 取均值并重命名为 VC
    control_data = {}
    if not vc_df.empty:
        control_data["VC"] = vc_df.mean(axis=1).values
        
    return control_data

def main():
    parser = argparse.ArgumentParser(description='处理 Thermal Overpress All 数据')
    parser.add_argument('--rawdata_dir', type=str, required=True,
                        help='原始数据目录，包含 ALL.xlsx')
    parser.add_argument('--plate_file', type=str, required=True,
                        help='映射文件路径，如 plate_gene_mapping.xlsx')
    parser.add_argument('--result_dir', type=str, required=True,
                        help='输出结果目录')
    parser.add_argument('--c1_label', type=str, default='stress',
                        choices=['stress', 'nonstress'],
                        help='Condition 1 标签')
    parser.add_argument('--c2_label', type=str, default='nonstress',
                        choices=['stress', 'nonstress'],
                        help='Condition 2 标签')
    
    args = parser.parse_args()
    
    rawdata_dir = Path(args.rawdata_dir)
    result_dir = Path(args.result_dir)
    result_dir.mkdir(parents=True, exist_ok=True)
    
    # 1. 加载映射
    position_to_gene = load_plate_gene_mapping(args.plate_file)
    print(f"已加载 {len(position_to_gene)} 个位置映射")
    
    # 2. 读取 ALL.xlsx
    all_excel_path = rawdata_dir / "ALL.xlsx"
    if not all_excel_path.exists():
        print(f"错误: 找不到文件 {all_excel_path}")
        return
    
    print(f"读取数据文件: {all_excel_path}")
    xls = pd.ExcelFile(all_excel_path)
    sheet_names = xls.sheet_names
    print(f"Sheets 列表: {sheet_names}")
    
    # 根据用户描述：
    # 第 1 个 sheet: 实验条件突变株 (Thermal_40)
    # 第 2 个 sheet: 对照条件突变株 (Control_30)
    # 第 3 个 sheet: VC 实验条件 (VC_40)
    # 第 4 个 sheet: VC 对照条件 (VC_30)
    
    # 处理 Sheet 1 (Experimental Mutant)
    df_mutant_c1 = pd.read_excel(xls, sheet_name=sheet_names[0])
    mutant_c1_data, time_points = process_mutant_sheet(df_mutant_c1, position_to_gene)
    
    # 处理 Sheet 2 (Control Mutant)
    df_mutant_c2 = pd.read_excel(xls, sheet_name=sheet_names[1])
    mutant_c2_data, _ = process_mutant_sheet(df_mutant_c2, position_to_gene)
    
    # 处理 Sheet 3 (VC Experimental)
    df_vc_c1 = pd.read_excel(xls, sheet_name=sheet_names[2])
    control_c1_data = process_vc_sheet(df_vc_c1)
    
    # 处理 Sheet 4 (VC Control)
    df_vc_c2 = pd.read_excel(xls, sheet_name=sheet_names[3])
    control_c2_data = process_vc_sheet(df_vc_c2)
    
    # 3. 保存 CSV 文件
    print("生成输出文件...")
    
    # 01. mutant_stress_OD.csv
    if mutant_c1_data:
        df = pd.DataFrame(mutant_c1_data, index=time_points)
        filename = f"01.mutant_{args.c1_label}_OD.csv"
        df.to_csv(result_dir / filename)
        print(f"已保存 {filename}")

    # 02. mutant_nonstress_OD.csv
    if mutant_c2_data:
        df = pd.DataFrame(mutant_c2_data, index=time_points)
        filename = f"02.mutant_{args.c2_label}_OD.csv"
        df.to_csv(result_dir / filename)
        print(f"已保存 {filename}")
        
    # 03. Ctrol_stress_OD.csv
    if control_c1_data:
        df = pd.DataFrame(control_c1_data, index=time_points)
        filename = f"03.Ctrol_{args.c1_label}_OD.csv"
        df.to_csv(result_dir / filename)
        print(f"已保存 {filename}")
        
    # 04. Ctrol_nonstress_OD.csv
    if control_c2_data:
        df = pd.DataFrame(control_c2_data, index=time_points)
        filename = f"04.Ctrol_{args.c2_label}_OD.csv"
        df.to_csv(result_dir / filename)
        print(f"已保存 {filename}")

if __name__ == "__main__":
    main()