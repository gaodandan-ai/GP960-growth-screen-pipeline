import pandas as pd
import matplotlib.pyplot as plt
import os
import argparse
from matplotlib.backends.backend_pdf import PdfPages

def plot_combined_curves(k_dir, o_dir, output_dir, time_max=None, mapping_file=None):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    # Load mapping if provided
    gene_map = {}
    if mapping_file and os.path.exists(mapping_file):
        try:
            df_map = pd.read_excel(mapping_file)
            # Ensure columns exist
            if 'old_locus_tag' in df_map.columns and 'Gene_Name' in df_map.columns:
                gene_map = pd.Series(df_map.Gene_Name.values, index=df_map.old_locus_tag).to_dict()
                print(f"Loaded {len(gene_map)} gene mappings.")
        except Exception as e:
            print(f"Error loading mapping file: {e}")
        
    conditions = ["stress", "nonstress"]
    
    for condition in conditions:
        print(f"Processing condition: {condition}")
        
        # Load mutant data
        k_file = os.path.join(k_dir, f"mutant_{condition}_OD_cleaned.csv")
        o_file = os.path.join(o_dir, f"mutant_{condition}_OD_cleaned.csv")
        
        if not os.path.exists(k_file) or not os.path.exists(o_file):
            print(f"Skipping {condition} because files are missing.")
            continue
            
        df_k = pd.read_csv(k_file)
        df_o = pd.read_csv(o_file)
        
        # Load control data
        ck_file = os.path.join(k_dir, f"Ctrol_{condition}_OD_cleaned.csv")
        co_file = os.path.join(o_dir, f"Ctrol_{condition}_OD_cleaned.csv")
        
        # Average controls for reference
        ctrl_k_avg = None
        if os.path.exists(ck_file):
            df_ck = pd.read_csv(ck_file)
            ctrl_k_avg = df_ck.iloc[:, 1:].mean(axis=1)
            
        ctrl_o_avg = None
        if os.path.exists(co_file):
            df_co = pd.read_csv(co_file)
            ctrl_o_avg = df_co.iloc[:, 1:].mean(axis=1)
        
        # Find common genes
        genes_k = df_k.columns[1:]
        genes_o = df_o.columns[1:]
        common_genes = sorted(list(set(genes_k).intersection(set(genes_o))))
        
        print(f"Found {len(common_genes)} common genes for {condition}")
        
        if not common_genes:
            continue

        # Plotting
        pdf_file = os.path.join(output_dir, f"combined_growth_curves_{condition}.pdf")
        with PdfPages(pdf_file) as pdf:
            # Group into grids of 3x3
            num_per_page = 9
            for i in range(0, len(common_genes), num_per_page):
                fig, axes = plt.subplots(3, 3, figsize=(15, 12))
                axes = axes.flatten()
                
                for j in range(num_per_page):
                    idx = i + j
                    if idx < len(common_genes):
                        gene = common_genes[idx]
                        ax = axes[j]
                        
                        time_k = df_k.iloc[:, 0] / 60.0 # Convert to hours
                        time_o = df_o.iloc[:, 0] / 60.0
                        
                        # Apply time_max if provided
                        if time_max is not None:
                            mask_k = time_k <= time_max
                            mask_o = time_o <= time_max
                        else:
                            mask_k = [True] * len(time_k)
                            mask_o = [True] * len(time_o)

                        ax.plot(time_k[mask_k], df_k[gene][mask_k], label='KO (K)', color='blue', linewidth=2)
                        ax.plot(time_o[mask_o], df_o[gene][mask_o], label='OE (O)', color='red', linewidth=2)
                        
                        # Plot controls
                        if ctrl_k_avg is not None:
                            ax.plot(time_k[mask_k], ctrl_k_avg[mask_k], label='Ctrl (K)', color='blue', linestyle='--', alpha=0.3)
                        if ctrl_o_avg is not None:
                            ax.plot(time_o[mask_o], ctrl_o_avg[mask_o], label='Ctrl (O)', color='red', linestyle='--', alpha=0.3)
                        
                        # Title with Gene Name if available
                        display_name = gene
                        if gene in gene_map and pd.notna(gene_map[gene]):
                            display_name = f"{gene} ({gene_map[gene]})"
                        
                        ax.set_title(f"Gene: {display_name}")
                        ax.set_xlabel("Time (h)")
                        ax.set_ylabel("OD")
                    else:
                        axes[j].axis('off')
                
                # Add a single legend for the entire figure at the bottom
                handles, labels = axes[0].get_legend_handles_labels()
                fig.legend(handles, labels, loc='lower center', ncol=4, fontsize='medium')
                
                # Adjust layout to make room for the legend at the bottom
                plt.tight_layout(rect=[0, 0.05, 1, 1])
                pdf.savefig(fig)
                plt.close(fig)
        
        print(f"Saved combined curves to {pdf_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Compare growth curves of common genes from two experiments.')
    parser.add_argument('--k_dir', type=str, required=True, help='Directory containing cleaned data for KO experiment')
    parser.add_argument('--o_dir', type=str, required=True, help='Directory containing cleaned data for OE experiment')
    parser.add_argument('--output_dir', type=str, required=True, help='Directory to save the results')
    parser.add_argument('--time_max', type=float, default=None, help='Maximum time to plot (in hours)')
    parser.add_argument('--mapping_file', type=str, default=None, help='Excel file containing gene mapping')
    
    args = parser.parse_args()
    
    plot_combined_curves(args.k_dir, args.o_dir, args.output_dir, args.time_max, args.mapping_file)
