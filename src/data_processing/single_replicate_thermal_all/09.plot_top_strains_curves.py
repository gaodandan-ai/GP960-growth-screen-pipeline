
#!/usr/bin/env python3
"""
12.plot_top_strains_curves.py
绘制各特征Top 5菌株的生长曲线。
每个特征的Top 5画在一起，即使有重复也重复画。
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from matplotlib.backends.backend_pdf import PdfPages
import warnings
warnings.filterwarnings('ignore')

# 设置图表样式
plt.rcParams['font.sans-serif'] = ['DejaVu Sans', 'Arial', 'SimHei']
plt.rcParams['axes.unicode_minus'] = False
sns.set_style("whitegrid")

def normalize_column_names(df, condition):
    """
    标准化列名，处理不同条件文件间的命名不一致问题
    """
    if condition in ['con_stress', 'con_nonstress']:
        new_df_data = {}
        for col in df.columns:
            if 'VC' in col:
                new_df_data[col] = df[col].values
            
            if '-' in col:
                gene_part = col.split('-', 1)[1]
                if ';' in gene_part:
                    genes = gene_part.split(';')
                    for gene in genes:
                        if gene.strip():
                            new_df_data[gene.strip()] = df[col].values
                else:
                    new_df_data[gene_part.strip()] = df[col].values
        new_df = pd.DataFrame(new_df_data, index=df.index)
        return new_df
    return df

def load_time_series_data(data_dir):
    """
    加载四种条件下的时间序列OD数据
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
            df = pd.read_csv(file_path, index_col=0)
            df = normalize_column_names(df, condition)
            # 时间已经是小时，不需要转换
            df.index = df.index
            data_dict[condition] = df
            print(f"Loaded {condition}: {df.shape}")
        else:
            print(f"Warning: {filename} not found")
    return data_dict

def plot_growth_curves_by_feature(data_dict, top5_file, output_dir, time_max, output_name='top_strains_curves_by_feature.pdf', gene_map={}):
    output_path = Path(output_dir) / output_name
    
    # 读取筛选结果
    df = pd.read_csv(top5_file)
    
    condition_labels = {
        'mut_stress': 'Mutant + Stress',
        'mut_nonstress': 'Mutant + Non-stress', 
        'con_stress': 'Control + Stress',
        'con_nonstress': 'Control + Non-stress'
    }
    
    condition_colors = {
        'mut_stress': '#d62728',      # Red
        'mut_nonstress': '#2ca02c',   # Green
        'con_stress': '#ff7f0e',      # Orange  
        'con_nonstress': '#1f77b4'    # Blue
    }
    
    # 按 Feature 分组
    features = df['Feature'].unique()
    
    with PdfPages(output_path) as pdf:
        for feature in features:
            feature_df = df[df['Feature'] == feature].sort_values('Rank')
            
            # 这一页画这个特征的 Top 5 (通常就是 5 个图)
            n_plots = len(feature_df)
            n_cols = 3
            n_rows = (n_plots + n_cols - 1) // n_cols
            
            fig, axes = plt.subplots(n_rows, n_cols, figsize=(18, 5 * n_rows), constrained_layout=True)
            
            # Handle single subplot case
            if n_rows == 1 and n_cols == 1:
                axes = np.array([axes])
            axes = axes.flatten() if isinstance(axes, np.ndarray) else [axes]
            
            # 添加大标题
            fig.suptitle(f"Top 5 Strains for Feature: {feature}", fontsize=16, fontweight='bold', y=1.02)
            
            for i, (_, row) in enumerate(feature_df.iterrows()):
                gene = row['Gene']
                rank = row['Rank']
                rf_val = row['RF_Value']
                
                ax = axes[i]
                
                # Plot lines
                curves_found = False
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
                                    label=label, linewidth=2, alpha=0.8)
                            curves_found = True
                
                # Title with Gene Name if available
                display_name = gene
                if gene in gene_map and pd.notna(gene_map[gene]):
                    display_name = f"{gene} ({gene_map[gene]})"
                
                # Title: Gene (Rank X)
                title = f"{display_name} (Rank {rank})\nRF = {rf_val:.2f}"
                ax.set_title(title, fontsize=12, fontweight='bold')
                ax.set_xlabel('Time (h)')
                ax.set_ylabel('OD600')
                ax.grid(True, alpha=0.3)
                ax.set_xlim(0, time_max)
                ax.set_ylim(0, None)
                
                if not curves_found:
                    ax.text(0.5, 0.5, "Data Not Found", ha='center', va='center')
                
            # Hide empty subplots
            for i in range(n_plots, len(axes)):
                axes[i].set_visible(False)
                
            # Legend (only once per page)
            handles, labels = axes[0].get_legend_handles_labels()
            if handles:
                fig.legend(handles, labels, loc='upper center', bbox_to_anchor=(0.5, 0.02), ncol=4, fontsize=12)
                
            pdf.savefig(fig, bbox_inches='tight')
            plt.close(fig)
            
    print(f"Saved curves to: {output_path}")

def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--input_dir', type=str, required=True, help='Path to cleaned data directory (02.cleaned_data)')
    parser.add_argument('--top5_file', type=str, required=True, help='Path to top5_strains_per_feature.csv')
    parser.add_argument('--output_dir', type=str, required=True)
    parser.add_argument('--time_max', type=float, default=25.0)
    parser.add_argument('--mapping_file', type=str, default=None, help='Excel file containing gene mapping')
    args = parser.parse_args()
    
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Load mapping
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
                print(f"Loaded {len(gene_map)} gene mappings")
        except Exception as e:
            print(f"Error loading mapping file: {e}")

    # 1. Load Data
    data_dict = load_time_series_data(args.input_dir)
    
    # 2. Plot
    plot_growth_curves_by_feature(data_dict, args.top5_file, output_dir, args.time_max, gene_map=gene_map)

if __name__ == "__main__":
    main()
