import pandas as pd
import matplotlib.pyplot as plt
import os
import argparse
from matplotlib.backends.backend_pdf import PdfPages

def plot_condition_comparison(k_dir, o_dir, output_dir, time_max=None, mapping_file=None):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        
    # Load mapping if provided
    gene_map = {}
    if mapping_file and os.path.exists(mapping_file):
        try:
            df_map = pd.read_excel(mapping_file)
            if 'old_locus_tag' in df_map.columns and 'Gene_Name' in df_map.columns:
                gene_map = pd.Series(df_map.Gene_Name.values, index=df_map.old_locus_tag).to_dict()
                print(f"Loaded {len(gene_map)} gene mappings.")
        except Exception as e:
            print(f"Error loading mapping file: {e}")

    # Load all necessary data
    data = {}
    for condition in ["stress", "nonstress"]:
        k_file = os.path.join(k_dir, f"mutant_{condition}_OD_cleaned.csv")
        o_file = os.path.join(o_dir, f"mutant_{condition}_OD_cleaned.csv")
        ck_file = os.path.join(k_dir, f"Ctrol_{condition}_OD_cleaned.csv")
        co_file = os.path.join(o_dir, f"Ctrol_{condition}_OD_cleaned.csv")
        
        if not all(os.path.exists(f) for f in [k_file, o_file]):
            print(f"Missing essential files for {condition}. Skipping.")
            continue
            
        df_k = pd.read_csv(k_file)
        df_o = pd.read_csv(o_file)
        
        ctrl_k_avg = None
        if os.path.exists(ck_file):
            df_ck = pd.read_csv(ck_file)
            ctrl_k_avg = df_ck.iloc[:, 1:].mean(axis=1)
            
        ctrl_o_avg = None
        if os.path.exists(co_file):
            df_co = pd.read_csv(co_file)
            ctrl_o_avg = df_co.iloc[:, 1:].mean(axis=1)
            
        data[condition] = {
            'df_k': df_k,
            'df_o': df_o,
            'ctrl_k': ctrl_k_avg,
            'ctrl_o': ctrl_o_avg,
            'time_k': df_k.iloc[:, 0] / 60.0,
            'time_o': df_o.iloc[:, 0] / 60.0
        }

    if "stress" not in data or "nonstress" not in data:
        print("Both stress and nonstress data are required.")
        return

    # Find common genes across all experiments
    genes_k = set(data['stress']['df_k'].columns[1:])
    genes_o = set(data['stress']['df_o'].columns[1:])
    genes_k_ns = set(data['nonstress']['df_k'].columns[1:])
    genes_o_ns = set(data['nonstress']['df_o'].columns[1:])
    
    common_genes = sorted(list(genes_k & genes_o & genes_k_ns & genes_o_ns))
    print(f"Found {len(common_genes)} common genes across conditions.")

    if not common_genes:
        return

    # Plotting
    pdf_file = os.path.join(output_dir, "gene_condition_comparison.pdf")
    with PdfPages(pdf_file) as pdf:
        # 4 genes per page, each has 2 plots (Stress/Non-stress) -> 4x2 grid
        num_genes_per_page = 4
        for i in range(0, len(common_genes), num_genes_per_page):
            fig, axes = plt.subplots(4, 2, figsize=(12, 16))
            
            for j in range(num_genes_per_page):
                idx = i + j
                if idx < len(common_genes):
                    gene = common_genes[idx]
                    
                    for col_idx, condition in enumerate(["nonstress", "stress"]):
                        ax = axes[j, col_idx]
                        d = data[condition]
                        
                        time_k = d['time_k']
                        time_o = d['time_o']
                        
                        if time_max is not None:
                            mask_k = time_k <= time_max
                            mask_o = time_o <= time_max
                        else:
                            mask_k = [True] * len(time_k)
                            mask_o = [True] * len(time_o)
                            
                        ax.plot(time_k[mask_k], d['df_k'][gene][mask_k], label='KO (K)', color='blue', linewidth=2)
                        ax.plot(time_o[mask_o], d['df_o'][gene][mask_o], label='OE (O)', color='red', linewidth=2)
                        
                        if d['ctrl_k'] is not None:
                            ax.plot(time_k[mask_k], d['ctrl_k'][mask_k], label='Ctrl (K)', color='blue', linestyle='--', alpha=0.3)
                        if d['ctrl_o'] is not None:
                            ax.plot(time_o[mask_o], d['ctrl_o'][mask_o], label='Ctrl (O)', color='red', linestyle='--', alpha=0.3)
                        
                        # Title with Gene Name if available
                        display_name = gene
                        if gene in gene_map and pd.notna(gene_map[gene]):
                            display_name = f"{gene} ({gene_map[gene]})"

                        ax.set_title(f"{display_name} ({condition})")
                        ax.set_xlabel("Time (h)")
                        ax.set_ylabel("OD")
                        ax.set_ylim(0, 22)
                        ax.set_yticks(range(0, 23, 5))
                else:
                    axes[j, 0].axis('off')
                    axes[j, 1].axis('off')
            
            # Shared legend at bottom
            handles, labels = axes[0, 0].get_legend_handles_labels()
            fig.legend(handles, labels, loc='lower center', ncol=4, fontsize='medium')
            
            plt.tight_layout(rect=[0, 0.03, 1, 0.97])
            pdf.savefig(fig)
            plt.close(fig)

    print(f"Saved condition comparisons to {pdf_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Compare stress vs nonstress conditions for common genes.')
    parser.add_argument('--k_dir', type=str, required=True, help='Directory containing cleaned data for KO experiment')
    parser.add_argument('--o_dir', type=str, required=True, help='Directory containing cleaned data for OE experiment')
    parser.add_argument('--output_dir', type=str, required=True, help='Directory to save the results')
    parser.add_argument('--time_max', type=float, default=None, help='Maximum time to plot (in hours)')
    parser.add_argument('--mapping_file', type=str, default=None, help='Excel file containing gene mapping')
    
    args = parser.parse_args()
    plot_condition_comparison(args.k_dir, args.o_dir, args.output_dir, args.time_max, args.mapping_file)
