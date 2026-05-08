#!/usr/bin/env python3
"""
04.plot_gene_od.py
绘制top10基因在四种条件下的时间-OD分布图

"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
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
        # Control条件文件的列名格式: ctrol1-Cgl0419;Cgl1715;...
        # 需要为每个基因创建单独的列
        new_df_data = {}
        
        for col in df.columns:
            if '-' in col and ';' in col:
                # 提取基因ID部分 (去掉ctrol1-前缀，分割分号分隔的基因)
                gene_part = col.split('-', 1)[1]  # 去掉ctrol1-前缀
                genes = gene_part.split(';')  # 分割基因
                
                # 为每个基因创建一个列，数据相同
                for gene in genes:
                    if gene.strip():  # 确保基因名不为空
                        new_df_data[gene.strip()] = df[col].values
            else:
                new_df_data[col] = df[col].values
        
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
        'mut_stress': '01.Goe_stress_OD_cleaned.csv',
        'mut_nonstress': '02.Goe_Nonstress_OD_cleaned.csv', 
        'con_stress': '03.Ctrol_stress_OD_cleaned.csv',
        'con_nonstress': '04.Ctrol_Nonstress_OD_cleaned.csv'
    }
    
    data_dict = {}
    
    for condition, filename in data_files.items():
        file_path = Path(data_dir) / filename
        if file_path.exists():
            # 所有文件现在都有时间列作为第一列
            df = pd.read_csv(file_path, index_col=0)
            # 标准化列名
            df = normalize_column_names(df, condition)
            # 将时间从分钟转换为小时
            df.index = df.index / 60.0
            
            data_dict[condition] = df
            print(f"已加载 {condition} 数据: {df.shape}, 时间范围: {df.index.min():.1f}-{df.index.max():.1f}小时")
        else:
            print(f"警告: 文件 {filename} 不存在")
    
    return data_dict

def get_top10_genes(result_dir, kind='stress'):
    """
    从指定CSV文件中获取top10基因列表
    
    Parameters:
    result_dir (str): 结果目录路径
    kind (str): 'stress' 或 'sensitive'
    
    Returns:
    list: top10基因名称列表
    """
    analysis_dir = Path(result_dir) / '04.optimized_analysis'
    
    if kind == 'stress':
        file_name = 'stress_tolerant_genes_optimized.csv'
    elif kind == 'sensitive':
        file_name = 'stress_sensitive_genes_optimized.csv'
    else:
        print(f"警告: 未知的基因类型 '{kind}'")
        return []

    file_path = analysis_dir / file_name
    
    if file_path.exists():
        df = pd.read_csv(file_path)
        top10_genes = df.head(10)['gene'].tolist()
        print(f"获取到 {kind} top 10 基因: {top10_genes}")
        return top10_genes
    else:
        print(f"警告: 文件 {file_path} 不存在")
        return []

def plot_gene_growth_curves(data_dict, gene_list, output_dir, title='Top 10 Genes Growth Curves', output_name='top10_genes_growth_curves.pdf'):
    """
    绘制指定基因在四种条件下的生长曲线
    
    Parameters:
    data_dict (dict): 包含四种条件数据的字典
    gene_list (list): 要绘制的基因列表
    output_dir (str): 输出目录路径
    title (str): 图标题
    output_name (str): 输出文件名（含扩展名）
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

    # 新增：条件样式和标记映射
    condition_styles = {
        'mut_stress':    {'linestyle': '-',  'marker': 'o'}, # 实线, 圆形
        'mut_nonstress': {'linestyle': '--', 'marker': 's'}, # 虚线, 方形
        'con_stress':    {'linestyle': '-.', 'marker': '^'}, # 点划线, 上三角
        'con_nonstress': {'linestyle': ':',  'marker': 'v'}  # 点线, 下三角
    }
    
    # 创建子图布局 (2行5列，显示10个基因)
    fig, axes = plt.subplots(2, 5, figsize=(20, 8))
    fig.suptitle(title, fontsize=16, fontweight='bold')
    
    axes = axes.flatten()
    
    for i, gene in enumerate(gene_list[:10]):  # 确保只绘制前10个基因
        ax = axes[i]
        
        # 为每种条件绘制生长曲线
        curves_plotted = 0
        for condition, label in condition_labels.items():
            if condition in data_dict and gene in data_dict[condition].columns:
                # 获取时间点和OD值
                time_points = data_dict[condition].index.values
                od_values = data_dict[condition][gene].values
                
                # 绘制曲线 (更新：增加样式和标记)
                ax.plot(time_points, od_values, 
                       color=condition_colors[condition], 
                       label=label, 
                       linewidth=0.5, 
                       alpha=0.8,
                       **condition_styles[condition],
                       markersize=1)
                curves_plotted += 1
                
            else:
                print(f"警告: 基因 {gene} 在条件 {condition} 中未找到")
        
        # 设置子图属性
        ax.set_title(f'{gene}', fontsize=16, fontweight='bold')
        ax.set_xlabel('Time (h)', fontsize=16,fontweight='bold') #Time
        ax.set_ylabel('OD600', fontsize=16,fontweight='bold') #OD600
        ax.grid(True, alpha=0.3)
        ax.set_xlim(0, 48)  # 设置为48小时
        ax.set_ylim(0, None)
        
        print(f"基因 {gene}: 绘制了 {curves_plotted} 条曲线")
    
    # 更新：在整个图的底部添加带样式的图例
    from matplotlib.lines import Line2D
    legend_elements = [
        Line2D([0], [0], 
               color=condition_colors[key], 
               label=condition_labels[key], 
               linewidth=1.5, 
               **condition_styles[key], 
               markersize=5)
        for key in condition_labels.keys()
    ]
    
    fig.legend(handles=legend_elements,
               loc='lower center', bbox_to_anchor=(0.5, -0.05), ncol=4, fontsize=20)
    
    # 调整布局，为底部图例留出空间
    plt.tight_layout()
    plt.subplots_adjust(bottom=0.1)
    
    # 保存为PDF
    output_path = Path(output_dir) / output_name
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"生长曲线图已保存至: {output_path}")
    
    plt.show()


def main():
    """主函数"""
    # 设置路径
    base_dir = Path('/data/zuoll/1.project/02.GOE/02_result/02.MSG')
    data_dir = base_dir / '02.cleaned_data'
    result_dir = base_dir
    output_dir = base_dir / '04.optimized_gene_growth_curves'
    
    # 创建输出目录
    output_dir.mkdir(exist_ok=True)
    
    print("=== 开始绘制Top10基因生长曲线 ===")
    
    # 1. 加载时间序列数据
    print("\n1. 加载时间序列数据...")
    data_dict = load_time_series_data(data_dir)
    
    if not data_dict:
        print("错误: 无法加载数据文件")
        return
    
    # 2. 获取stress与sensitive top10基因列表
    print("\n2. 获取top10基因列表...")
    stress_top10 = get_top10_genes(result_dir, kind='stress')
    sensitive_top10 = get_top10_genes(result_dir, kind='sensitive')
    
    if not stress_top10:
        print("错误: 无法获取stress top10基因列表")
    else:
        # 3. 绘制stress top10综合生长曲线图（PDF）
        print("\n3. 绘制stress top10综合生长曲线图...")
        plot_gene_growth_curves(
            data_dict, stress_top10, output_dir,
            title='', #Top 10 Stress-Tolerant Genes Growth Curves
            output_name='top10_stress_tolerant_genes_growth_curves_optimized.pdf'
        )

    if not sensitive_top10:
        print("错误: 无法获取sensitive top10基因列表")
    else:
        # 4. 绘制sensitive top10综合生长曲线图（PDF）
        print("\n4. 绘制sensitive top10综合生长曲线图...")
        plot_gene_growth_curves(
            data_dict, sensitive_top10, output_dir,
            title='', #Top 10 Stress-Sensitive Genes Growth Curves
            output_name='top10_stress_sensitive_genes_growth_curves_optimized.pdf'
        )
    
    # # 5. 绘制单个基因详细生长曲线图（保留原逻辑，仅对stress top10）
    # print("\n5. 绘制单个基因详细生长曲线图...")
    # plot_individual_gene_curves(data_dict, stress_top10, output_dir)
    
    # print("\n=== 生长曲线绘制完成 ===")
    # print(f"所有图片已保存至: {output_dir}")

    # 设置路径
    base_dir = Path('/data/zuoll/1.project/02.GOE/02_result/03.NaCl')
    data_dir = base_dir / '02.cleaned_data'
    result_dir = base_dir
    output_dir = base_dir / '04.optimized_gene_growth_curves'
    
    # 创建输出目录
    output_dir.mkdir(exist_ok=True)
    
    print("=== 开始绘制Top10基因生长曲线 ===")
    
    # 1. 加载时间序列数据
    print("\n1. 加载时间序列数据...")
    data_dict = load_time_series_data(data_dir)
    
    if not data_dict:
        print("错误: 无法加载数据文件")
        return
    
    # 2. 获取stress与sensitive top10基因列表
    print("\n2. 获取top10基因列表...")
    stress_top10 = get_top10_genes(result_dir, kind='stress')
    sensitive_top10 = get_top10_genes(result_dir, kind='sensitive')
    
    if not stress_top10:
        print("错误: 无法获取stress top10基因列表")
    else:
        # 3. 绘制stress top10综合生长曲线图（PDF）
        print("\n3. 绘制stress top10综合生长曲线图...")
        plot_gene_growth_curves(
            data_dict, stress_top10, output_dir,
            title='Top 10 Stress-Tolerant Genes Growth Curves (Optimized)',
            output_name='top10_stress_tolerant_genes_growth_curves_optimized.pdf'
        )

    if not sensitive_top10:
        print("错误: 无法获取sensitive top10基因列表")
    else:
        # 4. 绘制sensitive top10综合生长曲线图（PDF）
        print("\n4. 绘制sensitive top10综合生长曲线图...")
        plot_gene_growth_curves(
            data_dict, sensitive_top10, output_dir,
            title='Top 10 Stress-Sensitive Genes Growth Curves (Optimized)',
            output_name='top10_stress_sensitive_genes_growth_curves_optimized.pdf'
        )
    
    # # 5. 绘制单个基因详细生长曲线图（保留原逻辑，仅对stress top10）
    # print("\n5. 绘制单个基因详细生长曲线图...")
    # plot_individual_gene_curves(data_dict, stress_top10, output_dir)
    
    # print("\n=== 生长曲线绘制完成 ===")
    # print(f"所有图片已保存至: {output_dir}")

if __name__ == "__main__":
    main()