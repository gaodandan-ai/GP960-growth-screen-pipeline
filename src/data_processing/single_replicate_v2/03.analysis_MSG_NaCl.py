#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
NaCl NaCl胁迫耐受性分析
分析基因过表达实验中哪些基因在胁迫条件下表现出明显的耐受性。

分析指标：
- 生长性能指标 ：最大生长速率、最终生物量（OD最大值）、滞后期（OD达到初始值两倍的时间）、生长曲线下面积(AUC)
- 胁迫耐受性指标 ：相对适应性(Relative Fitness)、综合相对适应性评分
- 筛选标准 ：综合相对适应性 > 1.05, 最终OD相对适应性 > 1.05
- 非胁迫不耐受基因 ：综合相对适应性 < 0.95, 最终OD相对适应性 < 0.95

输出：
- 胁迫耐受性基因排名列表
- 非胁迫不耐受基因排名列表
- 生长指标统计表
- 可视化图表
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from scipy.optimize import curve_fit
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans, DBSCAN

# 设置matplotlib参数，避免中文显示问题
plt.rcParams['font.family'] = ['DejaVu Sans', 'sans-serif']
plt.rcParams['axes.unicode_minus'] = False

class StressToleranceAnalyzer:
    def __init__(self, data_dir):
        """初始化分析器"""
        self.data_dir = Path(data_dir)
        self.results = {}
        
    def load_data(self):
        """加载数据"""
        print("正在加载清理后的数据...")
        
        # 加载四个清理后的数据文件
        cleaned_data_dir = self.data_dir / "02.cleaned_data"
        self.goe_stress = pd.read_csv(cleaned_data_dir / "01.Goe_stress_OD_cleaned.csv", index_col=0)
        self.goe_nonstress = pd.read_csv(cleaned_data_dir / "02.Goe_Nonstress_OD_cleaned.csv", index_col=0)
        self.ctrl_stress = pd.read_csv(cleaned_data_dir / "03.Ctrol_stress_OD_cleaned.csv", index_col=0)
        self.ctrl_nonstress = pd.read_csv(cleaned_data_dir / "04.Ctrol_Nonstress_OD_cleaned.csv", index_col=0)
        
        # 获取时间点
        self.time_points = self.goe_stress.index.values
        
        # 获取基因列表（去除空格）
        self.genes = [col.strip() for col in self.goe_stress.columns]
        
        print(f"数据加载完成:")
        print(f"- 时间点: {len(self.time_points)} 个 ({self.time_points[0]:.1f} - {self.time_points[-1]:.1f} 分钟)")
        print(f"- 基因数: {len(self.genes)} 个")
        
    def logistic_growth(self, t, K, r, t0):
        """Logistic生长模型"""
        return K / (1 + np.exp(-r * (t - t0)))
    
    def calculate_growth_parameters(self, od_values, time_points):
        """计算生长参数"""
        try:
            # 移除异常值
            od_values = np.array(od_values)
            valid_idx = (od_values > 0) & (od_values < 20) & np.isfinite(od_values)
            
            if np.sum(valid_idx) < 10:  # 数据点太少
                return None
                
            od_clean = od_values[valid_idx]
            time_clean = time_points[valid_idx]
            
            # 拟合Logistic生长曲线
            try:
                # 初始参数估计
                K_init = np.max(od_clean)
                r_init = 0.01
                t0_init = np.median(time_clean)
                
                popt, _ = curve_fit(self.logistic_growth, time_clean, od_clean, 
                                  p0=[K_init, r_init, t0_init], 
                                  maxfev=1000, bounds=([0.1, 0.001, 0], [20, 1, 3000]))
                
                K, r, t0 = popt
                
            except:
                # 如果拟合失败，使用简单方法
                K = np.max(od_clean)
                r = 0.01
                t0 = time_clean[np.argmax(np.gradient(od_clean))]
            
            # 计算各项指标
            final_od = od_clean[-1]
            max_growth_rate = r * K / 4  # Logistic模型的最大生长速率
            
            # 滞后期：OD达到初始值2倍的时间
            lag_time = 0
            initial_od = od_clean[0]
            threshold = initial_od * 2
            lag_idx = np.where(od_clean >= threshold)[0]
            if len(lag_idx) > 0:
                lag_time = time_clean[lag_idx[0]]
            
            # 生长曲线下面积 (AUC)
            auc = np.trapz(od_clean, time_clean)
            
            return {
                'final_od': final_od,
                'max_growth_rate': max_growth_rate,
                'lag_time': lag_time,
                'auc': auc,
                'carrying_capacity': K,
                'growth_rate_constant': r,
                'inflection_point': t0
            }
            
        except Exception as e:
            print(f"计算生长参数时出错: {e}")
            return None
    
    def find_control_column(self, gene):
        """查找基因对应的对照组列名"""
        # 在对照组文件中查找包含该基因的列
        for col in self.ctrl_stress.columns:
            if gene in col:
                return col
        return None
    
    def analyze_gene_performance(self):
        """分析每个基因的生长性能（包含对照组数据计算相对适应性）"""
        print("正在分析基因生长性能（包含对照组数据）...")
        
        gene_results = []
        
        for i, gene in enumerate(self.genes):
            if i % 500 == 0:
                print(f"处理进度: {i+1}/{len(self.genes)}")
            
            try:
                # 获取该基因在过表达组的数据
                gene_col = gene + '  ' if gene + '  ' in self.goe_stress.columns else gene
                
                if gene_col not in self.goe_stress.columns:
                    continue
                
                # 获取过表达组数据
                mut_stress_data = self.goe_stress[gene_col].values
                mut_nonstress_data = self.goe_nonstress[gene_col].values
                
                # 查找对照组数据
                ctrl_col = self.find_control_column(gene)
                if ctrl_col is None:
                    continue
                
                con_stress_data = self.ctrl_stress[ctrl_col].values
                con_nonstress_data = self.ctrl_nonstress[ctrl_col].values
                
                # 计算生长参数
                mut_stress_params = self.calculate_growth_parameters(mut_stress_data, self.time_points)
                mut_nonstress_params = self.calculate_growth_parameters(mut_nonstress_data, self.time_points)
                con_stress_params = self.calculate_growth_parameters(con_stress_data, self.time_points)
                con_nonstress_params = self.calculate_growth_parameters(con_nonstress_data, self.time_points)
                
                if any(params is None for params in [mut_stress_params, mut_nonstress_params, 
                                                   con_stress_params, con_nonstress_params]):
                    continue
                
                # 计算相对适应性 (Relative Fitness)
                # RF = (OD_mut_stress / OD_mut_nonstress) / (OD_con_stress / OD_con_nonstress)
                mut_ratio = mut_stress_params['final_od'] / mut_nonstress_params['final_od'] if mut_nonstress_params['final_od'] > 0 else 0
                con_ratio = con_stress_params['final_od'] / con_nonstress_params['final_od'] if con_nonstress_params['final_od'] > 0 else 0
                relative_fitness = mut_ratio / con_ratio if con_ratio > 0 else 0
                
                # 计算基于AUC的相对适应性
                mut_auc_ratio = mut_stress_params['auc'] / mut_nonstress_params['auc'] if mut_nonstress_params['auc'] > 0 else 0
                con_auc_ratio = con_stress_params['auc'] / con_nonstress_params['auc'] if con_nonstress_params['auc'] > 0 else 0
                relative_fitness_auc = mut_auc_ratio / con_auc_ratio if con_auc_ratio > 0 else 0
                
                # 计算生长速率的相对适应性
                mut_rate_ratio = mut_stress_params['max_growth_rate'] / mut_nonstress_params['max_growth_rate'] if mut_nonstress_params['max_growth_rate'] > 0 else 0
                con_rate_ratio = con_stress_params['max_growth_rate'] / con_nonstress_params['max_growth_rate'] if con_nonstress_params['max_growth_rate'] > 0 else 0
                relative_fitness_rate = mut_rate_ratio / con_rate_ratio if con_rate_ratio > 0 else 0
                
                # 计算滞后期比例
                mut_lag_ratio = mut_stress_params['lag_time'] / mut_nonstress_params['lag_time'] if mut_nonstress_params['lag_time'] > 0 else 1
                con_lag_ratio = con_stress_params['lag_time'] / con_nonstress_params['lag_time'] if con_nonstress_params['lag_time'] > 0 else 1
                relative_lag_ratio = mut_lag_ratio / con_lag_ratio if con_lag_ratio > 0 else 1
                
                # 综合相对适应性评分 (权重: 最终OD 40%, AUC 40%, 生长速率 20%)
                comprehensive_rf = (relative_fitness * 0.4 + relative_fitness_auc * 0.4 + relative_fitness_rate * 0.2)
                
                # 胁迫耐受性评分 (考虑滞后期惩罚)
                lag_penalty = 1 / max(relative_lag_ratio, 0.1) if relative_lag_ratio > 1 else 1
                tolerance_score = comprehensive_rf * lag_penalty * 100
                
                gene_results.append({
                    'gene': gene,
                    'mut_stress_final_od': mut_stress_params['final_od'],
                    'mut_nonstress_final_od': mut_nonstress_params['final_od'],
                    'con_stress_final_od': con_stress_params['final_od'],
                    'con_nonstress_final_od': con_nonstress_params['final_od'],
                    'mut_stress_auc': mut_stress_params['auc'],
                    'mut_nonstress_auc': mut_nonstress_params['auc'],
                    'con_stress_auc': con_stress_params['auc'],
                    'con_nonstress_auc': con_nonstress_params['auc'],
                    'mut_stress_growth_rate': mut_stress_params['max_growth_rate'],
                    'mut_nonstress_growth_rate': mut_nonstress_params['max_growth_rate'],
                    'con_stress_growth_rate': con_stress_params['max_growth_rate'],
                    'con_nonstress_growth_rate': con_nonstress_params['max_growth_rate'],
                    'relative_fitness_od': relative_fitness,
                    'relative_fitness_auc': relative_fitness_auc,
                    'relative_fitness_rate': relative_fitness_rate,
                    'comprehensive_rf': comprehensive_rf,
                    'relative_lag_ratio': relative_lag_ratio,
                    'tolerance_score': tolerance_score,
                    'control_column': ctrl_col
                })
                
            except Exception as e:
                print(f"处理基因 {gene} 时出错: {e}")
                continue
        
        self.results_df = pd.DataFrame(gene_results)
        print(f"成功分析了 {len(self.results_df)} 个基因")
        
        # 新增：调用动态加权和优化CRF计算
        self.calculate_dynamic_weights_and_crf()
        
    def calculate_dynamic_weights_and_crf(self):
        """
        使用PCA动态计算各生长参数的权重，并生成优化的综合相对适应度(O-CRF)。
        """
        print("正在使用PCA计算动态权重和优化CRF...")
        
        # 1. 准备用于PCA的数据
        features = ['relative_fitness_od', 'relative_fitness_auc', 'relative_fitness_rate']
        pca_data = self.results_df[features].dropna()
        
        if pca_data.empty:
            print("警告: 没有足够的数据进行PCA分析。将使用固定权重。")
            self.results_df['optimized_tolerance_score'] = self.results_df['comprehensive_rf'] #optimized_tolerance_score
            self.results_df['optimized_tolerance_score'] = self.results_df['tolerance_score']
            self.dynamic_weights = [0.4, 0.4, 0.2]  # 保存固定权重
            return
            
        # 2. 数据标准化
        scaler = StandardScaler()
        scaled_data = scaler.fit_transform(pca_data)
        
        # 3. 执行PCA
        pca = PCA(n_components=1)
        pca.fit(scaled_data)
        
        # 4. 提取载荷作为动态权重
        # 载荷表示原始特征对主成分的贡献度
        loadings = pca.components_[0]
        
        # 将载荷转换为正权重（例如，取绝对值或平方后归一化）
        dynamic_weights = np.abs(loadings) / np.sum(np.abs(loadings))
        self.dynamic_weights = dynamic_weights  # 保存动态权重
        
        print(f"动态权重 (OD, AUC, Rate): {dynamic_weights[0]:.3f}, {dynamic_weights[1]:.3f}, {dynamic_weights[2]:.3f}")
        
        # 5. 使用动态权重计算Optimized CRF (O-CRF)
        self.results_df['optimized_crf'] = (
            self.results_df['relative_fitness_od'] * dynamic_weights[0] +
            self.results_df['relative_fitness_auc'] * dynamic_weights[1] +
            self.results_df['relative_fitness_rate'] * dynamic_weights[2]
        ) #optimized_crf
        
        # 6. 计算Optimized Tolerance Score (O-TS)
        lag_penalty = 1 / self.results_df['relative_lag_ratio'].clip(lower=0.1)
        lag_penalty[self.results_df['relative_lag_ratio'] <= 1] = 1
        self.results_df['optimized_tolerance_score'] = self.results_df['optimized_crf'] * lag_penalty * 100 #optimized_crf
        
        print("动态加权和优化CRF计算完成。")

    def filter_genes_with_kmeans(self, n_clusters=3):
        """
        使用K-Means在二维指标(O-TS 和 relative_fitness_od)上实现自适应阈值筛选。
        """
        print(f"正在使用K-Means在二维(O-TS, OD)上进行自适应阈值筛选 (n_clusters={n_clusters})...")
        
        # 1. 准备用于聚类的数据
        features = ['optimized_tolerance_score', 'relative_fitness_od'] #relative_fitness_od
        cluster_data = self.results_df[features].dropna()
        
        if len(cluster_data) < n_clusters:
            print("警告: 数据点不足以进行K-Means聚类。")
            return pd.DataFrame(), pd.DataFrame(), {}

        # 2. 数据标准化（二维聚类必须进行标准化）
        scaler = StandardScaler()
        scaled_features = scaler.fit_transform(cluster_data)

        # 3. 执行K-Means聚类
        kmeans = KMeans(n_clusters=n_clusters, random_state=42, n_init=10)
        self.results_df['cluster'] = kmeans.fit_predict(scaled_features)
        
        # 4. 识别“高耐受”和“敏感”簇
        # 策略：计算每个簇中心在原始空间的坐标，根据 O-CRF 和 OD 的综合表现排序
        # 这里简单取每个簇的 O-CRF 和 OD 的均值之和
        cluster_stats = self.results_df.groupby('cluster')[features].mean()
        cluster_performance = cluster_stats.sum(axis=1)
        
        tolerant_cluster_label = cluster_performance.idxmax()
        sensitive_cluster_label = cluster_performance.idxmin()
        
        # 5. 筛选出相应簇的基因
        tolerant_genes = self.results_df[self.results_df['cluster'] == tolerant_cluster_label].copy()
        sensitive_genes = self.results_df[self.results_df['cluster'] == sensitive_cluster_label].copy()
        
        # 6. 生成动态阈值 (基于两个指标的最小值/最大值)
        adaptive_thresholds = {
            'tolerant_threshold_ts': tolerant_genes['optimized_tolerance_score'].min(),
            'tolerant_threshold_od': tolerant_genes['relative_fitness_od'].min(), #relative_fitness_od
            'sensitive_threshold_ts': sensitive_genes['optimized_tolerance_score'].max(),
            'sensitive_threshold_od': sensitive_genes['relative_fitness_od'].max() #relative_fitness_od
        }
        
        print(f"自适应阈值确定 (二维):")
        print(f"- 胁迫耐受性边界: O-TS > {adaptive_thresholds['tolerant_threshold_ts']:.3f}, OD > {adaptive_thresholds['tolerant_threshold_od']:.3f}")
        print(f"- 胁迫敏感性边界: O-TS < {adaptive_thresholds['sensitive_threshold_ts']:.3f}, OD < {adaptive_thresholds['sensitive_threshold_od']:.3f}")
        
        # 7. 排序
        tolerant_genes = tolerant_genes.sort_values('optimized_tolerance_score', ascending=False)
        sensitive_genes = sensitive_genes.sort_values('optimized_tolerance_score', ascending=True)
        
        print(f"筛选出 {len(tolerant_genes)} 个胁迫耐受性基因 (K-Means 2D)")
        print(f"筛选出 {len(sensitive_genes)} 个胁迫敏感性基因 (K-Means 2D)")
        
        return tolerant_genes, sensitive_genes, adaptive_thresholds

    
    def generate_summary_statistics(self):
        """生成汇总统计"""
        print("生成汇总统计...")
        
        summary = {
            'total_genes': len(self.results_df),
            'relative_fitness_od_stats': self.results_df['relative_fitness_od'].describe(),
            'relative_fitness_auc_stats': self.results_df['relative_fitness_auc'].describe(),
            'comprehensive_rf_stats': self.results_df['comprehensive_rf'].describe(),
            'tolerance_score_stats': self.results_df['tolerance_score'].describe()
        }

        # Add new stats if they exist
        if 'optimized_tolerance_score' in self.results_df.columns:
            summary['optimized_tolerance_score_stats'] = self.results_df['optimized_tolerance_score'].describe()
        if 'optimized_tolerance_score' in self.results_df.columns:
            summary['optimized_tolerance_score_stats'] = self.results_df['optimized_tolerance_score'].describe()
        
        return summary
    
    def plot_pca_weights(self, output_dir):
        """绘制PCA动态权重的条形图"""
        if not hasattr(self, 'dynamic_weights'):
            return

        print("绘制PCA动态权重图...")
        features = ['RF_OD', 'RF_AUC', 'RF_Rate']
        weights = self.dynamic_weights

        fig, ax = plt.subplots(figsize=(8, 6))
        ax.bar(features, weights, color=['#3498db', '#2ecc71', '#e74c3c'])
        ax.set_ylabel('Dynamic Weight')
        ax.set_title('PCA-Derived Dynamic Weights for CRF Calculation')
        ax.set_ylim(0, 1.0)

        for i, w in enumerate(weights):
            ax.text(i, w + 0.02, f'{w:.3f}', ha='center', va='bottom')

        fig.tight_layout()
        fig.savefig(output_dir / "pca_dynamic_weights.pdf", dpi=300)
        plt.close(fig)

    def plot_tolerance_distribution(self, output_dir):
        """绘制耐受性指标分布图及二维聚类散点图"""
        print("绘制耐受性指标分布图及二维聚类散点图...")
        
        # 1) 优化后的综合相对适应性分布 (按K-Means簇着色)
        fig, ax = plt.subplots(1, 1, figsize=(10, 7))
        sns.histplot(data=self.results_df, x='optimized_tolerance_score', hue='cluster', 
                     palette='viridis', multiple="stack", ax=ax, bins=100, legend=True)
        
        # 提取自适应阈值
        if hasattr(self, 'adaptive_thresholds') and self.adaptive_thresholds:
            t_threshold_ts = self.adaptive_thresholds.get('tolerant_threshold_ts')
            s_threshold_ts = self.adaptive_thresholds.get('sensitive_threshold_ts')
            if t_threshold_ts:
                ax.axvline(t_threshold_ts, color='red', linestyle='--', label=f'Tolerant O-TS ({t_threshold_ts:.2f})')
            if s_threshold_ts:
                ax.axvline(s_threshold_ts, color='blue', linestyle='--', label=f'Sensitive O-TS ({s_threshold_ts:.2f})')
        
        ax.set_xlabel('Optimized Comprehensive Relative Fitness (O-TS)')
        ax.set_ylabel('Number of Genes')
        ax.set_title('K-Means Clustering Distribution on Optimized TS')
        ax.legend()
        fig.tight_layout()
        fig.savefig(output_dir / "optimized_tolerance_score_distribution_kmeans.pdf", dpi=300, bbox_inches='tight')
        plt.close(fig)

        # 2) 二维聚类散点图: O-CRF vs OD
        fig, ax = plt.subplots(figsize=(10, 8))
        sns.scatterplot(data=self.results_df, x='optimized_tolerance_score', y='relative_fitness_od', 
                        hue='cluster', palette='viridis', style='cluster', s=60, alpha=0.7, ax=ax)
        
        # 绘制自适应阈值线
        if hasattr(self, 'adaptive_thresholds') and self.adaptive_thresholds:
            t_ts = self.adaptive_thresholds.get('tolerant_threshold_ts')
            t_od = self.adaptive_thresholds.get('tolerant_threshold_od')
            s_ts = self.adaptive_thresholds.get('sensitive_threshold_ts')
            s_od = self.adaptive_thresholds.get('sensitive_threshold_od')
            
            # 耐受性边界
            if t_ts and t_od:
                ax.axvline(t_ts, color='red', linestyle='--', alpha=0.6)
                ax.axhline(t_od, color='red', linestyle='--', alpha=0.6, label='Tolerant Boundary')
            
            # 敏感性边界
            if s_ts and s_od:
                ax.axvline(s_ts, color='blue', linestyle='--', alpha=0.6)
                ax.axhline(s_od, color='blue', linestyle='--', alpha=0.6, label='Sensitive Boundary')

        ax.set_xlabel('Optimized Tolerance Score (O-TS)')
        ax.set_ylabel('Relative Fitness (OD)')
        ax.set_title('2D K-Means Adaptive Clustering (O-TS vs OD)')
        ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        fig.tight_layout()
        fig.savefig(output_dir / "adaptive_clustering_scatter.pdf", dpi=300, bbox_inches='tight')
        plt.close(fig)
        
    def save_results(self, output_dir):
        """保存分析结果"""
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        
        print(f"保存结果到 {output_dir}")
        
        # 使用K-Means方法筛选基因
        tolerant_genes, sensitive_genes, adaptive_thresholds = self.filter_genes_with_kmeans(n_clusters=3)
        self.adaptive_thresholds = adaptive_thresholds # 保存阈值以便绘图
        
        # 保存完整结果
        self.results_df.to_csv(output_dir / "all_genes_analysis_optimized.csv", index=False)
        
        # 保存筛选出的基因
        tolerant_genes.to_csv(output_dir / "stress_tolerant_genes_optimized.csv", index=False)
        sensitive_genes.to_csv(output_dir / "stress_sensitive_genes_optimized.csv", index=False)
        
        # 保存汇总统计
        summary = self.generate_summary_statistics()
        with open(output_dir / "summary_statistics_optimized.txt", 'w', encoding='utf-8') as f:
            f.write("胁迫耐受性分析汇总报告 (PCA动态加权 + 1D-KMeans自适应阈值)\n")
            f.write("=" * 70 + "\n\n")
            f.write(f"总基因数: {summary['total_genes']}\n\n")

            if hasattr(self, 'dynamic_weights'):
                f.write(f"PCA动态权重 (OD, AUC, Rate): {self.dynamic_weights[0]:.3f}, {self.dynamic_weights[1]:.3f}, {self.dynamic_weights[2]:.3f}\n\n")
            
            if adaptive_thresholds:
                f.write(f"2D-KMeans自适应阈值 (O-TS, OD):\n")
                f.write(f"- 胁迫耐受性边界: O-TS > {adaptive_thresholds.get('tolerant_threshold_ts', 'N/A'):.3f}, OD > {adaptive_thresholds.get('tolerant_threshold_od', 'N/A'):.3f}\n")
                f.write(f"- 胁迫敏感性边界: O-TS < {adaptive_thresholds.get('sensitive_threshold_ts', 'N/A'):.3f}, OD < {adaptive_thresholds.get('sensitive_threshold_od', 'N/A'):.3f}\n\n")
            
            f.write(f"胁迫耐受性基因数 (K-Means): {len(tolerant_genes)} ({len(tolerant_genes)/summary['total_genes']*100:.2f}%)\n")
            f.write(f"胁迫敏感性基因数 (K-Means): {len(sensitive_genes)} ({len(sensitive_genes)/summary['total_genes']*100:.2f}%)\n\n")
            
            if 'optimized_tolerance_score_stats' in summary:
                f.write("优化综合相对适应性 (O-CRF) 统计:\n")
                f.write(str(summary['optimized_tolerance_score_stats']) + "\n\n")
            
            if 'optimized_tolerance_score_stats' in summary:
                f.write("优化胁迫耐受性评分 (O-TS) 统计:\n")
                f.write(str(summary['optimized_tolerance_score_stats']) + "\n\n")

            f.write("原始综合相对适应性统计:\n")
            f.write(str(summary['comprehensive_rf_stats']) + "\n\n")

            f.write("相对适应性 (最终OD) 统计:\n")
            f.write(str(summary['relative_fitness_od_stats']) + "\n\n")
        
        # 生成可视化图表
        self.plot_pca_weights(output_dir)
        self.plot_tolerance_distribution(output_dir)
        
        print("结果保存完成!")
        return tolerant_genes

def main():
    """主函数"""
    # --- MSG 分析 ---
    print("\n" + "="*80)
    print("开始处理 MSG 胁迫数据集")
    print("="*80)
    msg_data_dir = Path("/data/zuoll/1.project/02.GOE/02_result/02.MSG")
    msg_output_dir = Path("/data/zuoll/1.project/02.GOE/02_result/02.MSG/04.optimized_analysis")
    
    analyzer_msg = StressToleranceAnalyzer(msg_data_dir)
    analyzer_msg.load_data()
    analyzer_msg.analyze_gene_performance()
    tolerant_genes_msg = analyzer_msg.save_results(msg_output_dir)
    
    print("\n" + "-"*60)
    print("MSG 胁迫 - 前10个最佳耐受性基因 (优化后):")
    print("-"*60)
    if not tolerant_genes_msg.empty:
        top_genes_msg = tolerant_genes_msg.head(10)
        for i, (_, gene_data) in enumerate(top_genes_msg.iterrows(), 1):
            print(f"{i:2d}. {gene_data['gene']}")
            print(f"    优化综合评分 (O-TS): {gene_data['optimized_tolerance_score']:.2f}")
            print(f"    优化综合适应性 (O-CRF): {gene_data['optimized_tolerance_score']:.3f}")
            print(f"    相对适应性(OD): {gene_data['relative_fitness_od']:.3f}")
    else:
        print("未找到符合条件的胁迫耐受性基因。")
    print(f"\n详细结果已保存到: {msg_output_dir}")

    # --- NaCl 分析 ---
    print("\n" + "="*80)
    print("开始处理 NaCl 胁迫数据集")
    print("="*80)
    nacl_data_dir = Path("/data/zuoll/1.project/02.GOE/02_result/03.NaCl")
    nacl_output_dir = Path("/data/zuoll/1.project/02.GOE/02_result/03.NaCl/04.optimized_analysis")
    
    analyzer_nacl = StressToleranceAnalyzer(nacl_data_dir)
    analyzer_nacl.load_data()
    analyzer_nacl.analyze_gene_performance()
    tolerant_genes_nacl = analyzer_nacl.save_results(nacl_output_dir)
    
    print("\n" + "-"*60)
    print("NaCl 胁迫 - 前10个最佳耐受性基因 (优化后):")
    print("-"*60)
    if not tolerant_genes_nacl.empty:
        top_genes_nacl = tolerant_genes_nacl.head(10)
        for i, (_, gene_data) in enumerate(top_genes_nacl.iterrows(), 1):
            print(f"{i:2d}. {gene_data['gene']}")
            print(f"    优化综合评分 (O-TS): {gene_data['optimized_tolerance_score']:.2f}")
            print(f"    优化综合适应性 (O-CRF): {gene_data['optimized_tolerance_score']:.3f}")
            print(f"    相对适应性(OD): {gene_data['relative_fitness_od']:.3f}")
    else:
        print("未找到符合条件的胁迫耐受性基因。")
    print(f"\n详细结果已保存到: {nacl_output_dir}")

if __name__ == "__main__":
    main()