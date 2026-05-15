# GP960 生长筛选分析流水线 (Growth Screen Pipeline)

[English Version](./README.md)

这是一个基于 Python 的分析流水线，专门用于处理由 GP960 系统生成的高通量微生物生长曲线筛选数据。

该流水线能够将原始的板式 Excel 文件转换为标准化的 OD 表格，清理低质量测量数据，计算胁迫相关的适应性指标（Fitness），提取生长特征，并生成可视化摘要，用于筛选耐受性菌株或基因。

## 主要功能

- **多设计支持**：支持单重复（Single-replicate）和三重复（Triplicate）实验设计的原始 `.xlsx` 文件处理。
- **数据清理**：使用可配置的 OD 阈值清理异常或低生长曲线。
- **适应性指标**：计算相对适应性 (RF)、基于 AUC 的 RF 以及综合胁迫耐受性评分。
- **深度特征提取**：提取 50 多个生长特征，包括最大 OD、AUC、生长速率、倍增时间以及 Gompertz 模型参数。
- **可视化**：生成生长曲线图、相关性热图、特征散点图、Top 菌株图和比较点图。

## 工作原理

### 1. 数据标准化
流水线将 GP960 原始 Excel 文件（时间点为列，孔位为行）解析为“Tidy”格式：
- **时间**：作为索引（行）。
- **孔位/基因**：作为列。
- 来自多个板的数据将被合并为每个实验条件（如胁迫 vs 非胁迫）的单个 CSV 文件。

### 2. 数据清理与筛选
为了确保分析的鲁棒性，流水线会自动过滤低质量的生长曲线：
- **基于阈值**：在非胁迫条件下未达到最小 OD（例如 OD > 1.0）的菌株将被标记。
- **异常检测**：具有过度噪声或非生物性跳跃的曲线将从高级适应性排名中剔除。

### 3. 适应性计算 (Relative Fitness, RF)
胁迫耐受性的核心衡量指标是**相对适应性 (RF)**。它将突变体的胁迫生长与其自身的非胁迫生长进行归一化，并与对照组（如空载对照 VC）进行比较：

$$RF = \frac{OD_{mutant, stress} / OD_{mutant, nonstress}}{OD_{control, stress} / OD_{control, nonstress}}$$

除了终点 OD，流水线还会计算：
- **AUC RF**：基于曲线下面积。
- **Rate RF**：基于最大比生长速率。
- **综合 RF (Comprehensive RF)**：加权评分（通常为 40% OD + 40% AUC + 20% Rate）。
- **耐受性评分 (Tolerance Score)**：综合 RF 经**迟滞期惩罚 (Lag Phase Penalty)** ($1/LagRatio$) 修正后的得分。

### 4. 生长特征提取
利用 `scipy.optimize` 和 `numpy`，流水线可以提取 50 多个特征：
- **Logistic/Gompertz 拟合**：提取环境容纳量 ($K$)、生长速率 ($r$) 和拐点。
- **Mu 计算**：通过滑动窗口对数线性回归计算最大比生长速率。
- **阶段检测**：自动识别迟滞期、对数期和稳定期。

## 仓库结构

```text
GP960-growth-screen-pipeline/
|-- src/
|   |-- data_processing/
|   |   |-- single_replicate/    # 标准单重复工作流
|   |   `-- triple_replicate/    # 标准三重复工作流
|   `-- pipline/                 # 用于各种实验的 Bash 批处理脚本
|-- config/                      # 实验设计映射和配置
`-- data/                        # 本地数据（Git 已忽略）
    |-- raw/                     # 原始 GP960 Excel 文件
    `-- results/                 # 生成的分析结果
```

## 安装

```bash
git clone https://github.com/gaodandan-ai/GP960-growth-screen-pipeline.git
cd GP960-growth-screen-pipeline
conda env create -f environmental.yml
conda activate gp960_analysis
```

## 工作流示例 (单重复)

### 1. 预处理 (Preprocess)
```bash
python src/data_processing/single_replicate/01.pp_rawdata.py \
    --rawdata_dir data/raw/highmethanol \
    --plate_file data/raw/plate-gene-mapping.xlsx \
    --result_dir data/results/01.ppraw_data \
    --c1_label stress --c2_label nonstress
```

### 2. 数据清理 (Clean)
```bash
python src/data_processing/single_replicate/02.data_cleaning.py \
    --input_dir data/results/01.ppraw_data \
    --output_dir data/results/02.cleaned_data \
    --mapping_file data/raw/plate-gene-mapping.xlsx \
    --threshold 1.0
```

### 3. 适应性分析 (Analyze Fitness)
```bash
python src/data_processing/single_replicate/03.calculate_fitness.py \
    --input_dir data/results/02.cleaned_data \
    --output_dir data/results/03.fitness_analysis \
    --rf_threshold_stress 1.05
```

## 输出结构

- `01.ppraw_data/`: 合并后的原始 OD 表格。
- `03.fitness_analysis/`: `stress_tolerant_genes.csv` (排名列表) 及分布图。
- `05.full_features/`: `all_features_raw_50plus.csv` (主特征表)。
- `09.plot_top5_curves_by_feature/`: 表现最佳候选菌株的生长曲线 PDF。

## 注意事项

- **时间单位**：请确保各脚本间时间单位的一致性（分钟 vs 小时）。
- **通用对照组**：部分高级脚本支持使用单个 "VC" 列作为所有板子的通用对照。
- **目录兼容性**：目录名 `src/pipline/` 为保持系统兼容性而沿用。

## 许可证

尚未指定许可证。
