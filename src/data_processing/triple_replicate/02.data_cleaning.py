#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
数据清理脚本
清理AmGlu胁迫耐受性分析中的异常数据

清理规则：
1. 删除plate-gene-mapping.xlsx中available=0基因
2. 对于所有文件: 先计算三个生物学重复的平均OD,作为该基因的OD值
3. 删除平均OD或最大OD小于某一个阈值的基因，如OD<1
4. 对于01.02文件: 删除异常基因所在列
5. 对于03.04文件: 只计算平均值，不删除基因列
"""

import pandas as pd
import numpy as np
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

class DataCleaner:
    def __init__(self, data_dir, output_dir,mapping_file):
        """初始化数据清理器"""
        self.data_dir = Path(data_dir)
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        # 数据文件路径
        self.file_paths = {
            f'mutant_nonstress': self.data_dir / f"mutant_nonstress_OD.csv",
            f'mutant_stress': self.data_dir / f"mutant_stress_OD.csv", 
            f'ctrl_nonstress': self.data_dir / f"Ctrol_nonstress_OD.csv",
            f'ctrl_stress': self.data_dir / f"Ctrol_stress_OD.csv"
        }
        
        # 基因映射文件
        self.mapping_file = mapping_file
        
        self.genes_to_remove = set()
        
    def load_gene_mapping(self):
        """加载基因映射文件，识别available=False的基因"""
        print("加载基因映射文件...")
        
        try:
            # 读取Excel文件
            mapping_df = pd.read_excel(self.mapping_file)
            print(f"基因映射文件列名: {mapping_df.columns.tolist()}")
            print(f"基因映射文件形状: {mapping_df.shape}")
            
            # 查看前几行数据
            print("\n前5行数据:")
            print(mapping_df.head())
            
            # 查找available列（可能有不同的列名，包括拼写错误的avalible）
            available_col = None
            for col in mapping_df.columns:
                if 'available' in col.lower() or 'avail' in col.lower() or 'avalible' in col.lower():
                    available_col = col
                    break
            
            if available_col is None:
                print("未找到available列，查看所有列:")
                for col in mapping_df.columns:
                    print(f"  {col}: {mapping_df[col].dtype}")
                return set()
            
            print(f"找到available列: {available_col}")
            
            # 查找基因ID列（优先选择gene列）
            gene_col = None
            for col in mapping_df.columns:
                if 'gene' in col.lower():
                    gene_col = col
                    break
            
            if gene_col is None:
                # 如果没有找到gene列，再查找包含id的列
                for col in mapping_df.columns:
                    if 'id' in col.lower():
                        gene_col = col
                        break
            
            if gene_col is None:
                # 假设第一列是基因ID
                gene_col = mapping_df.columns[0]
            
            print(f"使用基因ID列: {gene_col}")
            
            # 筛选available=0的基因（不可用基因）
            unavailable_mask = mapping_df[available_col] == 0
            unavailable_genes = mapping_df[unavailable_mask][gene_col].tolist()
            
            # 确保基因ID是字符串类型，并去重
            unavailable_genes = list(set([str(gene) for gene in unavailable_genes if pd.notna(gene)]))
            
            print(f"找到 {len(unavailable_genes)} 个不可用基因")
            if len(unavailable_genes) > 0:
                print(f"前10个不可用基因: {sorted(unavailable_genes)[:10]}")
            
            return set(unavailable_genes)
            
        except Exception as e:
            print(f"读取基因映射文件时出错: {e}")
            return set()
    
    def load_od_data(self):
        """加载OD数据文件"""
        print("加载OD数据文件...")
        
        self.od_data = {}
        
        for name, file_path in self.file_paths.items():
            if file_path.exists():
                df = pd.read_csv(file_path, index_col=0)
                self.od_data[name] = df
                print(f"{name}: {df.shape}")
            else:
                print(f"文件不存在: {file_path}")
    
    def calculate_replicate_averages(self):
        """计算所有文件中三个生物学重复的平均值"""
        print("计算三个生物学重复的平均值...")
        
        # 处理所有文件
        for file_key in self.od_data.keys():
            print(f"\n处理 {file_key}:")
            df = self.od_data[file_key]
            
            # 获取所有列名
            columns = df.columns.tolist()
            
            if file_key in ['mutant_stress', 'mutant_nonstress']:
                # Mutant文件：标准格式 xxxx-1, xxxx-2, xxxx-3
                gene_groups = {}
                for col in columns:
                    if '-' in col: #if col.startswith('Cgl') and '-' in col:
                        # 提取基因名部分（去掉-1、-2、-3后缀）
                        gene_name = col.rsplit('-', 1)[0]  # 从右边分割，取第一部分
                        if gene_name not in gene_groups:
                            gene_groups[gene_name] = []
                        gene_groups[gene_name].append(col)
                
                print(f"  发现 {len(gene_groups)} 个基因")
                
                # 创建新的DataFrame存储平均值
                averaged_data = pd.DataFrame(index=df.index)
                
                # 对每个基因计算三个重复的平均值
                for gene_name, replicate_cols in gene_groups.items():
                    if len(replicate_cols) == 3:  # 确保有三个重复
                        # 计算三个重复的平均值
                        avg_values = df[replicate_cols].mean(axis=1)
                        # 新列名不带编号，只保留基因名
                        averaged_data[gene_name] = avg_values
                    else:
                        print(f"  警告: 基因 {gene_name} 的重复数不是3个: {replicate_cols}")
                        # 如果不是3个重复，仍然计算平均值
                        avg_values = df[replicate_cols].mean(axis=1)
                        averaged_data[gene_name] = avg_values
            
            elif file_key in ['ctrl_stress', 'ctrl_nonstress']:
                # Control文件：特殊格式 ctrol1-genes;genes;..., ctrol2-genes;genes;..., ctrol3-genes;genes;...
                # 按重复分组列名
                replicate_groups = {'ctrol1-': [], 'ctrol2-': [], 'ctrol3-': []}
                
                for col in columns:
                    for prefix in replicate_groups.keys():
                        if col.startswith(prefix):
                            replicate_groups[prefix].append(col)
                            break
                
                print(f"  发现重复组: ctrol1-({len(replicate_groups['ctrol1-'])}), ctrol2-({len(replicate_groups['ctrol2-'])}), ctrol3-({len(replicate_groups['ctrol3-'])})")
                
                # 创建新的DataFrame存储平均值
                averaged_data = pd.DataFrame(index=df.index)
                
                # 找到最大列数
                max_cols = max(len(cols) for cols in replicate_groups.values())
                
                for i in range(max_cols):
                    # 获取对应位置的三个重复列
                    replicate_cols = []
                    for prefix in ['ctrol1-', 'ctrol2-', 'ctrol3-']:
                        if i < len(replicate_groups[prefix]):
                            replicate_cols.append(replicate_groups[prefix][i])
                    
                    if len(replicate_cols) == 3:  # 确保有三个重复
                        # 计算三个重复的平均值
                        avg_values = df[replicate_cols].mean(axis=1)
                        # 使用第一个重复的列名（去掉ctrol1-前缀）
                        col_name = replicate_groups['ctrol1-'][i].replace('ctrol1-', '')
                        averaged_data[col_name] = avg_values
                    else:
                        print(f"  警告: 位置 {i} 只有 {len(replicate_cols)} 个重复")
            
            # 更新od_data中的数据
            self.od_data[file_key] = averaged_data
            
            print(f"  原始列数: {len(columns)}")
            print(f"  平均后列数: {len(averaged_data.columns)}")
            print(f"  处理完成")
                
    def identify_low_od_genes(self, threshold=1.0):
        """识别在任一条件下平均OD或最大OD小于threshold的基因"""
        print(f"识别平均OD或最大OD < {threshold}的基因...")
        
        low_od_genes = set()
        
        for condition_name, data in self.od_data.items():
            print(f"\n分析 {condition_name}:")
            
            # 计算每个基因的平均OD和最大OD
            avg_od = data.mean(axis=0)
            max_od = data.max(axis=0)
            
            # 找出平均OD或最大OD小于阈值的基因
            low_avg_genes = avg_od[avg_od < threshold].index.tolist()
            low_max_genes = max_od[max_od < threshold].index.tolist()
            
            condition_low_genes = set(low_avg_genes + low_max_genes)
            low_od_genes.update(condition_low_genes)
            
            print(f"  平均OD < {threshold}的基因: {len(low_avg_genes)}")
            print(f"  最大OD < {threshold}的基因: {len(low_max_genes)}")
            print(f"  该条件下异常基因总数: {len(condition_low_genes)}")
            
            if len(condition_low_genes) > 0:
                print(f"  前5个异常基因: {list(condition_low_genes)[:5]}")
        
        print(f"\n所有条件下异常基因总数: {len(low_od_genes)}")
        return low_od_genes
    
    def clean_goe_files(self, genes_to_remove):
        """清理01.02文件（GOE文件），直接删除整列"""
        print("清理GOE文件（01.02）...")
        
        goe_files = ['mutant_stress', 'mutant_nonstress']
        
        for file_key in goe_files:
            if file_key in self.od_data:
                original_data = self.od_data[file_key].copy()
                
                # 找出需要删除的列
                cols_to_remove = [col for col in original_data.columns if col in genes_to_remove]
                
                # 删除列
                cleaned_data = original_data.drop(columns=cols_to_remove)
                
                print(f"{file_key}:")
                print(f"  原始基因数: {len(original_data.columns)}")
                print(f"  删除基因数: {len(cols_to_remove)}")
                print(f"  剩余基因数: {len(cleaned_data.columns)}")
                
                # 保存清理后的数据
                output_file = self.output_dir / f"{self.file_paths[file_key].name.replace('.csv', '_cleaned.csv')}"
                cleaned_data.to_csv(output_file)
                print(f"  保存到: {output_file}")
    
    def clean_control_files(self):
        """处理03,04的control文件，保存三个生物学重复的平均值"""
        print("处理Control文件...")
        
        control_files = ['ctrl_stress', 'ctrl_nonstress']
        
        for file_key in control_files:
            if file_key in self.od_data:
                print(f"\n处理 {file_key}:")
                
                # 获取已经计算好平均值的数据
                averaged_data = self.od_data[file_key]
                
                # 生成输出文件路径
                output_file = self.output_dir / f"{self.file_paths[file_key].name.replace('.csv', '_cleaned.csv')}"
                
                # 保存平均后的数据
                averaged_data.to_csv(output_file)
                
                print(f"  列数: {len(averaged_data.columns)}")
                print(f"  保存到: {output_file}")
                print(f"  处理完成")
            else:
                print(f"警告: {file_key} 数据未加载")
    
    def generate_cleaning_report(self, unavailable_genes, low_od_genes, total_genes_to_remove):
        """生成清理报告"""
        report_file = self.output_dir / "data_cleaning_report.txt"
        
        with open(report_file, 'w', encoding='utf-8') as f:
            f.write("数据清理报告\n")
            f.write("=" * 50 + "\n\n")
            
            f.write("清理规则:\n")
            f.write("1. 删除plate-gene-mapping.xlsx中available=False的基因\n")
            f.write("2. 对所有文件计算三个生物学重复的平均值\n")
            f.write("3. 删除在任一条件下平均OD或最大OD小于1的基因\n")
            f.write("4. GOE文件删除异常基因列，Control文件保留所有基因\n\n")
            
            f.write("清理结果:\n")
            f.write(f"不可用基因数量: {len(unavailable_genes)}\n")
            f.write(f"低OD基因数量: {len(low_od_genes)}\n")
            f.write(f"总删除基因数量: {len(total_genes_to_remove)}\n\n")
            
            if len(unavailable_genes) > 0:
                f.write("不可用基因列表:\n")
                for gene in sorted(unavailable_genes):
                    f.write(f"  {gene}\n")
                f.write("\n")
            
            if len(low_od_genes) > 0:
                f.write("低OD基因列表（前50个）:\n")
                for gene in sorted(list(low_od_genes)[:50]):
                    f.write(f"  {gene}\n")
                f.write("\n")
            
            f.write("文件处理情况:\n")
            for name, file_path in self.file_paths.items():
                if name in self.od_data:
                    final_cols = len(self.od_data[name].columns)
                    if name in ['goe_stress', 'goe_nonstress']:
                        f.write(f"{file_path.name}: 最终基因数 {final_cols}\n")
                    else:
                        f.write(f"{file_path.name}: 最终列数 {final_cols}\n")
        
        print(f"清理报告保存到: {report_file}")
    
    def run_cleaning(self, threshold):
        """执行完整的数据清理流程"""
        print("开始数据清理流程...")
        
        # 1. 加载基因映射，识别不可用基因
        unavailable_genes = self.load_gene_mapping()
        
        # 2. 加载OD数据
        self.load_od_data()
        
        # 3. 计算所有文件中三个生物学重复的平均值
        self.calculate_replicate_averages()
        
        # 4. 识别低OD基因
        low_od_genes = self.identify_low_od_genes(threshold)
        
        # 5. 合并所有需要删除的基因
        total_genes_to_remove = unavailable_genes | low_od_genes
        
        print(f"\n总计需要删除的基因数量: {len(total_genes_to_remove)}")
        
        # 6. 清理GOE文件（01.02）
        self.clean_goe_files(total_genes_to_remove)
        
        # 7. 清理Control文件（03.04）- 计算三个重复的平均值
        self.clean_control_files()
        
        # 8. 生成清理报告
        self.generate_cleaning_report(unavailable_genes, low_od_genes, total_genes_to_remove)
        
        print("\n数据清理完成!")
        return total_genes_to_remove

def main():
    """主函数"""
    import argparse
    parser = argparse.ArgumentParser(description='处理高通量生长曲线数据')
    parser.add_argument('--input_dir', type=str, required=True,
                        help='输入数据目录(01.ppraw_data)')
    parser.add_argument('--output_dir', type=str, required=True,
                        help='输出数据目录(02.cleaned_data)')
    parser.add_argument('--mapping_file', type=str, required=True,
                        help='基因映射文件路径')
    parser.add_argument('--threshold', type=float, required=True,
                        help='处理的OD值阈值，如: 1')
    args = parser.parse_args()
    
    # 创建数据清理器
    cleaner = DataCleaner(args.input_dir, args.output_dir, args.mapping_file)
    
    try:
        # 执行清理
        removed_genes = cleaner.run_cleaning(args.threshold)
        
        print(f"\n清理完成!")
        print(f"删除了 {len(removed_genes)} 个异常基因（仅从GOE文件中删除）")
        print(f"清理后的数据保存在: {args.output_dir}")
        
    except Exception as e:
        print(f"数据清理过程中出现错误: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main()