#!/usr/bin/env python3
"""
04.plot_gene_od.py
绘制all基因在四种条件下的时间-OD分布图

"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from matplotlib.backends.backend_pdf import PdfPages
import warnings
warnings.filterwarnings('ignore')

# 设置中文字体和图表样式
plt.rcParams['font.sans-serif'] = ['DejaVu Sans', 'SimHei', 'Arial Unicode MS']
plt.rcParams['axes.unicode_minus'] = False
sns.set_style("whitegrid")

def normalize_column_names(df, condition):
    """
    标准化列名，处理不同条件文件间的命名不一致问题
    
    Parameters:
    df (DataFrame): 原始数据框
    condition (str): 条件名称
    
    Returns:
    DataFrame: 标准化列名后的数据框
    """
    if condition in ['con_stress', 'con_nonstress']:
        # Control条件文件的列名格式: ctrol1-Cgl0419;Cgl1715;... 或 ctrol-VC
        new_df_data = {}
        vc_data = None
        
        for col in df.columns:
            if 'VC' in col:
                vc_data = df[col].values
                new_df_data[col] = vc_data
            
            if '-' in col:
                # 提取基因ID部分
                gene_part = col.split('-', 1)[1]
                if ';' in gene_part:
                    genes = gene_part.split(';')
                    for gene in genes:
                        if gene.strip():
                            new_df_data[gene.strip()] = df[col].values
                else:
                    # 单个基因或 VC
                    new_df_data[gene_part.strip()] = df[col].values
        
        # 如果有 VC 数据，但某些基因没有对应列，可以在这里处理
        # 但通常绘图时会检查，所以这里保证 VC 列存在即可
        
        # 创建新的DataFrame
        new_df = pd.DataFrame(new_df_data, index=df.index)
        print(f"  标准化 {condition} 列名: 从 {len(df.columns)} 列扩展到 {len(new_df.columns)} 列")
        return new_df
    
    return df

def load_time_series_data(data_dir):
    """
    加载四种条件下的时间序列OD数据
    
    Parameters:
    data_dir (str): 数据目录路径
    
    Returns:
    dict: 包含四种条件数据的字典
    """
    data_files = {
        'mut_stress': '01.mutant_stress_OD_cleaned.csv',
        'mut_nonstress': '02.mutant_nonstress_OD_cleaned.csv', 
        'con_stress': '03.Ctrol_stress_OD_cleaned.csv',
        'con_nonstress': '04.Ctrol_nonstress_OD_cleaned.csv'
    }
    
    data_dict = {}
    
    for condition, filename in data_files.items():
        file_path = Path(data_dir) / filename
        if file_path.exists():
            # 所有文件现在都有时间列作为第一列
            df = pd.read_csv(file_path, index_col=0)
            # 标准化列名
            df = normalize_column_names(df, condition)
            # 时间已经是小时，不需要转换
            df.index = df.index
            
            data_dict[condition] = df
            print(f"已加载 {condition} 数据: {df.shape}, 时间范围: {df.index.min():.1f}-{df.index.max():.1f}小时")
        else:
            print(f"警告: 文件 {filename} 不存在")
    
    return data_dict

def get_all_genes(result_dir, kind='stress'):
    """
    从指定CSV文件中获取all基因列表
    
    Parameters:
    result_dir (str): 结果目录路径
    kind (str): 'stress' 或 'nonstress'
    
    Returns:
    list: all基因名称列表
    """
    analysis_dir = Path(result_dir)
    file_path = analysis_dir / ('stress_tolerant_genes.csv' if kind == 'stress' else 'nonstress_tolerant_genes.csv')
    
    if file_path.exists():
        df = pd.read_csv(file_path)
        all_genes = df['gene'].tolist()
        print(f"获取到{kind} all基因: {all_genes}")
        return all_genes
    else:
        print(f"警告: 文件 {file_path} 不存在")
        return []

def plot_gene_growth_curves(data_dict, gene_list, output_dir, time_max, title='Genes Growth Curves', output_name='genes_growth_curves.pdf', chunk_size=20, n_cols_max=5, dpi=150, gene_map={}):
    """
    绘制指定基因在四种条件下的生长曲线（动态子图布局）
    
    Parameters:
    data_dict (dict): 包含四种条件数据的字典
    gene_list (list): 要绘制的基因列表
    output_dir (str): 输出目录路径
    title (str): 图标题
    output_name (str): 输出文件名（含扩展名）
    gene_map (dict): 基因ID到名称的映射
    """
    # 条件标签和颜色映射
    condition_labels = {
        'mut_stress': 'Mutant + Stress',
        'mut_nonstress': 'Mutant + Non-stress', 
        'con_stress': 'Control + Stress',
        'con_nonstress': 'Control + Non-stress'
    }
    
    condition_colors = {
        'mut_stress': '#d62728',      # 红色
        'mut_nonstress': '#2ca02c',   # 绿色
        'con_stress': '#ff7f0e',      # 橙色  
        'con_nonstress': '#1f77b4'    # 蓝色
    }
    
    n_genes = len(gene_list)
    if n_genes == 0:
        print("警告: 没有基因数据可绘制")
        return
    output_path = Path(output_dir) / output_name
    with PdfPages(output_path) as pdf:
        for start in range(0, n_genes, chunk_size):
            subset = gene_list[start:start + chunk_size]
            n_cols = min(n_cols_max, len(subset))
            n_rows = (len(subset) + n_cols - 1) // n_cols
            fig_height = max(6, 3.5 * n_rows)
            fig, axes = plt.subplots(n_rows, n_cols, figsize=(20, fig_height), constrained_layout=True)
            fig.suptitle(title, fontsize=16, fontweight='bold')
            if n_rows == 1 and n_cols == 1:
                axes = [axes]
            elif n_rows == 1:
                axes = axes
            else:
                axes = axes.flatten()
            for i, gene in enumerate(subset):
                ax = axes[i]
                curves_plotted = 0
                for condition, label in condition_labels.items():
                    if condition in data_dict:
                        col_name = None
                        if not condition.startswith('con_'):
                            # 突变株组找基因名
                            if gene in data_dict[condition].columns:
                                col_name = gene
                        else:
                            # 控制组统一找 VC
                            for col in data_dict[condition].columns:
                                if 'VC' in col:
                                    col_name = col
                                    break
                        
                        if col_name:
                            time_points = data_dict[condition].index.values
                            od_values = data_dict[condition][col_name].values
                            ax.plot(time_points, od_values,
                                    color=condition_colors[condition],
                                    label=label,
                                    linewidth=2,
                                    alpha=0.8)
                            curves_plotted += 1
                
                # 设置标题，包含基因名（如果有）
                display_name = gene
                if gene in gene_map and pd.notna(gene_map[gene]):
                    display_name = f"{gene} ({gene_map[gene]})"
                
                ax.set_title(f'{display_name}', fontsize=12, fontweight='bold')
                ax.set_xlabel('Time (h)', fontsize=10)
                ax.set_ylabel('OD600', fontsize=10)
                ax.grid(True, alpha=0.3)
                ax.set_xlim(0, time_max)
                ax.set_ylim(0, None)
                print(f"基因 {gene}: 绘制了 {curves_plotted} 条曲线")
            for i in range(len(subset), len(axes)):
                axes[i].set_visible(False)
            fig.legend(labels=['Mutant + Stress', 'Mutant + Non-stress', 'Control + Stress', 'Control + Non-stress'],
                       loc='lower center', bbox_to_anchor=(0.5, -0.02), ncol=4, fontsize=12)
            pdf.savefig(fig, dpi=dpi, bbox_inches='tight')
            plt.close(fig)
    print(f"生长曲线图已保存至: {output_path} (共{n_genes}个基因，分页输出)")


def main():
    """主函数"""
    import argparse
    parser = argparse.ArgumentParser(description='处理高通量生长曲线数据')
    parser.add_argument('--input_dir', type=str, required=True,
                        help='输入数据目录')
    parser.add_argument('--output_dir', type=str, required=True,
                        help='输出数据目录(04.gene_growth_curves)')
    parser.add_argument('--analysis_dir', type=str, required=True,
                        help='分析结果目录(03.stress_tolerance_analysis)')
    parser.add_argument('--time_max', type=float, default=25.0,
                        help='绘制的时间范围(小时)，默认25')
    parser.add_argument('--mapping_file', type=str, default=None,
                        help='基因映射文件路径 (Excel)')
    args = parser.parse_args()

    # 设置路径
    data_dir = Path(args.input_dir)
    result_dir = Path(args.analysis_dir)
    output_dir = Path(args.output_dir)
    
    # 加载映射
    gene_map = {}
    if args.mapping_file and Path(args.mapping_file).exists():
        try:
            df_map = pd.read_excel(args.mapping_file)
            tag_col = None
            name_col = None
            for col in df_map.columns:
                if col.lower() in ['old_locus_tag', 'gene', 'id']: tag_col = col
                if col.lower() in ['gene_name', 'name']: name_col = col
            if tag_col and name_col:
                gene_map = pd.Series(df_map[name_col].values, index=df_map[tag_col]).to_dict()
                print(f"加载了 {len(gene_map)} 个基因映射")
        except Exception as e:
            print(f"加载映射文件失败: {e}")

    # 创建输出目录
    output_dir.mkdir(parents=True, exist_ok=True)
    
    print("=== 开始绘制all基因生长曲线 ===")
    
    # 1. 加载时间序列数据
    print("\n1. 加载时间序列数据...")
    data_dict = load_time_series_data(data_dir)
    
    if not data_dict:
        print("错误: 无法加载数据文件")
        return
    
    # 2. 获取stress与nonstress all基因列表
    print("\n2. 获取all基因列表...")
    stress_all = get_all_genes(result_dir, kind='stress')
    nonstress_all = get_all_genes(result_dir, kind='nonstress')
    
    if not stress_all:
        print("错误: 无法获取stress all基因列表")
        return
    if not nonstress_all:
        print("错误: 无法获取nonstress all基因列表")
        return
    
    # 3. 绘制stress all综合生长曲线图（PDF）
    print("\n3. 绘制stress all综合生长曲线图...")
    plot_gene_growth_curves(
        data_dict, stress_all, output_dir, args.time_max,
        title='ALL Stress-Tolerant Genes Growth Curves',
        output_name='all_stress_tolerant_genes_growth_curves.pdf',
        gene_map=gene_map
    )
    
    # 4. 绘制nonstress all综合生长曲线图（PDF）
    print("\n4. 绘制nonstress all综合生长曲线图...")
    plot_gene_growth_curves(
        data_dict, nonstress_all, output_dir, args.time_max,
        title='ALL Non-stress Tolerant Genes Growth Curves',
        output_name='all_nonstress_tolerant_genes_growth_curves.pdf',
        gene_map=gene_map
    )
    
    print("\n=== 生长曲线绘制完成 ===")
    print(f"所有图片已保存至: {output_dir}")

if __name__ == "__main__":
    main()
