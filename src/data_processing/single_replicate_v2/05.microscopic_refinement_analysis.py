
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from scipy.signal import savgol_filter
import os

class MicroscopicFeatureAnalyzer:
    def __init__(self, base_dir, analysis_name):
        self.base_dir = Path(base_dir)
        self.analysis_name = analysis_name
        self.data_dir = self.base_dir / "02.cleaned_data"
        self.macro_analysis_dir = self.base_dir / "04.optimized_analysis"
        self.output_dir = self.base_dir / "05.microscopic_refinement"
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
    def load_data(self):
        """加载宏观筛选结果和原始OD数据"""
        print(f"正在加载 {self.analysis_name} 的数据...")
        
        # 1. 加载一级筛选（K-Means）结果
        macro_file = self.macro_analysis_dir / "stress_tolerant_genes_optimized.csv"
        if not macro_file.exists():
            print(f"错误: 找不到一级筛选结果文件 {macro_file}")
            return False
        self.candidate_genes_df = pd.read_csv(macro_file)
        
        # 2. 加载过表达组胁迫条件下的原始OD数据
        od_file = self.data_dir / "01.Goe_stress_OD_cleaned.csv"
        if not od_file.exists():
            print(f"错误: 找不到OD数据文件 {od_file}")
            return False
        self.od_data = pd.read_csv(od_file, index_col=0)
        self.time_points = self.od_data.index.values / 60.0  # 转换为小时
        
        print(f"成功加载 {len(self.candidate_genes_df)} 个候选基因。")
        return True

    def extract_micro_dynamics(self, gene_name):
        """二级微观形态校准：提取高阶导数特征"""
        if gene_name not in self.od_data.columns:
            # 尝试处理带空格的情况
            gene_name_alt = gene_name + '  '
            if gene_name_alt in self.od_data.columns:
                od_values = self.od_data[gene_name_alt].values
            else:
                return None
        else:
            od_values = self.od_data[gene_name].values
            
        # 1. 数据平滑
        try:
            smooth_od = savgol_filter(od_values, window_length=11, polyorder=3)
            
            # 2. 计算一阶导数（生长速率 velocity）
            dt = np.mean(np.diff(self.time_points))
            velocity = np.gradient(smooth_od, dt)
            
            # 3. 计算二阶导数（生长加速度 acceleration）
            acceleration = np.gradient(velocity, dt)
            
            # 4. 提取核心微观特征
            max_v = np.max(velocity)  # 最大瞬时速率
            max_a = np.max(acceleration)  # 最大瞬时加速度（爆发力）
            
            # 生长稳定性评估：计算对数生长期后的速率衰减平滑度
            v_max_idx = np.argmax(velocity)
            deceleration_phase = velocity[v_max_idx:]
            if len(deceleration_phase) > 5:
                # 衰减抖动度：二阶导数在衰减期的波动标准差
                stability_noise = np.std(np.gradient(deceleration_phase, dt))
            else:
                stability_noise = 999
                
            return {
                'gene': gene_name,
                'max_velocity': max_v,
                'max_acceleration': max_a,
                'stability_noise': stability_noise,
                'smooth_od': smooth_od,
                'velocity_curve': velocity,
                'acceleration_curve': acceleration
            }
        except Exception as e:
            print(f"分析基因 {gene_name} 的微观动力学时出错: {e}")
            return None

    def perform_refinement(self):
        """执行多级递进式筛选"""
        print("执行二级微观形态校准筛选...")
        
        refined_results = []
        all_dynamics = {}
        
        for _, row in self.candidate_genes_df.iterrows():
            gene = row['gene']
            dynamics = self.extract_micro_dynamics(gene)
            
            if dynamics:
                # 独立验证算子设计：
                # 1. 加速度（爆发力）需处于前列
                # 2. 稳定性噪声（生长衰退期的波动）需低于平均水平
                # 这里我们先计算所有候选基因的平均指标
                all_dynamics[gene] = dynamics
                
        # 计算微观指标阈值（基于候选群体的统计分布）
        avg_noise = np.median([d['stability_noise'] for d in all_dynamics.values() if d['stability_noise'] < 100])
        avg_accel = np.median([d['max_acceleration'] for d in all_dynamics.values()])
        
        for gene, d in all_dynamics.items():
            # 筛选标准：爆发力不低于中位数的 80%，且生长过程抖动（噪声）不高于中位数的 1.5 倍
            # 这能有效剔除那些虽然长得多但生长极不稳定的“假阳性”
            is_micro_stable = (d['max_acceleration'] >= avg_accel * 0.8) and (d['stability_noise'] <= avg_noise * 1.5)
            
            # 将微观特征合并入结果
            macro_row = self.candidate_genes_df[self.candidate_genes_df['gene'] == gene].iloc[0].to_dict()
            macro_row.update({
                'micro_max_acceleration': d['max_acceleration'],
                'micro_stability_noise': d['stability_noise'],
                'is_micro_stable': 'Yes' if is_micro_stable else 'No'
            })
            refined_results.append(macro_row)
            
        self.refined_df = pd.DataFrame(refined_results)
        # 按照宏观和微观综合排序：(标准化后的 O-TS + 标准化后的加速度)
        self.refined_df['refinement_rank_score'] = (
            (self.refined_df['optimized_tolerance_score'] / self.refined_df['optimized_tolerance_score'].max()) + 
            (self.refined_df['micro_max_acceleration'] / self.refined_df['micro_max_acceleration'].max())
        )
        self.refined_df = self.refined_df.sort_values('refinement_rank_score', ascending=False)
        
        # 保存结果
        self.refined_df.to_csv(self.output_dir / "refined_tolerant_genes_micro_calibration.csv", index=False)
        print(f"精细化筛选完成。筛选出稳定耐受基因: {len(self.refined_df[self.refined_df['is_micro_stable'] == 'Yes'])}")
        
        return all_dynamics

    def plot_micro_dynamics(self, all_dynamics, top_n=3):
        """可视化微观动力学特征对比 (使用不同线型以适配黑白打印)"""
        print(f"筛选并绘制前 {top_n} 个稳健耐受基因的微观动力学图表...")
        
        # 首先过滤出 is_micro_stable 为 'Yes' 的基因，然后再取前 top_n
        stable_genes_df = self.refined_df[self.refined_df['is_micro_stable'] == 'Yes']
        
        if stable_genes_df.empty:
            print("警告: 没有基因通过稳健性校准，跳过绘图。")
            return
            
        top_genes = stable_genes_df.head(top_n)['gene'].tolist()
        
        fig, axes = plt.subplots(3, 1, figsize=(12, 15), sharex=True)
        
        # 定义不同的线型和标记以区分基因（适配黑白）
        linestyles = ['-', '--', ':', '-.', (0, (3, 5, 1, 5))]
        markers = ['o', 's', '^', 'D', 'v']
        colors = sns.color_palette("colorblind", top_n) # 使用色盲友好调色板，黑白转换效果更好
        
        for i, gene in enumerate(top_genes):
            d = all_dynamics[gene]
            ls = linestyles[i % len(linestyles)]
            mk = markers[i % len(markers)]
            
            # 1. 生长曲线 (OD)
            axes[0].plot(self.time_points, d['smooth_od'], label=f"{gene}", 
                        color=colors[i], linestyle=ls, linewidth=2)
            # 2. 生长速率 (1st Derivative)
            axes[1].plot(self.time_points, d['velocity_curve'], label=f"{gene}", 
                        color=colors[i], linestyle=ls, linewidth=2)
            # 3. 生长加速度 (2nd Derivative)
            axes[2].plot(self.time_points, d['acceleration_curve'], label=f"{gene}", 
                        color=colors[i], linestyle=ls, linewidth=1.5, alpha=0.8)
            
            # 每隔一段距离画一个标记点，增强黑白区分度
            markevery = max(1, len(self.time_points) // 10)
            axes[0].plot(self.time_points, d['smooth_od'], color=colors[i], 
                        marker=mk, markersize=6, markevery=markevery, linestyle='None')
            
        axes[0].set_ylabel("OD", fontsize=16, fontweight='bold')
        axes[0].set_title("", fontsize=14)
        axes[1].set_ylabel("dOD/dt", fontsize=16, fontweight='bold')
        axes[1].set_title("", fontsize=14)
        axes[2].set_ylabel("d²OD/dt²", fontsize=16, fontweight='bold')
        axes[2].set_title("", fontsize=14)
        axes[2].set_xlabel("Time (h)", fontsize=16, fontweight='bold')
        
        for ax in axes:
            ax.legend(loc='upper left', bbox_to_anchor=(1, 1), fontsize=16)
            ax.grid(True, alpha=0.3)
            # 设置坐标轴刻度数值的字体大小
            ax.tick_params(axis='both', labelsize=14)
            
        plt.tight_layout()
        plt.savefig(self.output_dir / "micro_refinement_dynamics_comparison.pdf", dpi=300, bbox_inches='tight')
        plt.close()

def main():
    # 对 MSG 数据集进行精细化筛选
    msg_analyzer = MicroscopicFeatureAnalyzer(
        base_dir="/data/zuoll/1.project/02.GOE/02_result/02.MSG",
        analysis_name="MSG_Stress"
    )
    
    if msg_analyzer.load_data():
        dynamics = msg_analyzer.perform_refinement()
        msg_analyzer.plot_micro_dynamics(dynamics)
        
    # 对 NaCl 数据集进行精细化筛选
    nacl_analyzer = MicroscopicFeatureAnalyzer(
        base_dir="/data/zuoll/1.project/02.GOE/02_result/03.NaCl",
        analysis_name="NaCl_Stress"
    )
    
    if nacl_analyzer.load_data():
        dynamics = nacl_analyzer.perform_refinement()
        nacl_analyzer.plot_micro_dynamics(dynamics)

if __name__ == "__main__":
    main()
