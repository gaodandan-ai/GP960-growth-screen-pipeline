#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
计算相对适应性矩阵 (RF Matrix)

RF = (Mut_Stress / Mut_NonStress) / (Con_Stress / Con_NonStress)

对于"越低越好"的特征（如滞后期、衰减幅度），取倒数 1/RF，
使得 RF 统一变为"越高 = 该基因在胁迫下表现越好"。

输入: all_features_raw_50plus.csv (由 05 脚本生成)
      selected_features_list.csv  (含 Feature 列)
输出: rf_matrix_final.csv
"""

import pandas as pd
import numpy as np
from pathlib import Path
import argparse

# ------------------------------------------------------------------
# 特征名须与 05.calculate_full_features.py 的输出列名精确一致。
# "越低越好"的特征：RF 计算后取倒数，使高 RF 始终代表更好表现。
# ------------------------------------------------------------------
LOWER_IS_BETTER = {
    # 时间类：越短越好
    'lag_phase',
    't_detect',
    't_lag_end',
    't10', 't50', 't90',
    'time_max_OD',        # 05 实际输出名 (已从 "time_max_OD / t_peak" 改为此)
    't_peak_from_start',
    'Lag_cost_index',
    # 衰减类：越小越好
    'decline_fraction',
    'max_decline_rate',   # 通常为负数，1/RF 逻辑见下方说明
    # 拟合质量：越小越好
    'fit_RMSE',
    'AIC',
    'BIC',
    # 噪声：越小越好
    'post_plateau_noise',
}

# max_decline_rate 通常为负值处理说明：
#   Mut=-0.1 (衰减慢,好), Con=-0.5 (衰减快,差) → ratio=0.2 → 1/ratio=5 (高=好) ✓
#   Mut=-0.5 (衰减快,差), Con=-0.1 (衰减慢,好) → ratio=5   → 1/ratio=0.2 (低=差) ✓

EPSILON = 1e-9   # 防止除零，比旧代码的 1e-6 更严格


def _safe_rf(ms: pd.Series, mn: pd.Series,
             cs: pd.Series, cn: pd.Series,
             feature: str) -> pd.Series:
    """
    向量化计算单个特征的 RF。
    RF = (ms/mn) / (cs/cn)
    若 feature 在 LOWER_IS_BETTER 中则取 1/RF。
    """
    # 将接近 0 的分母替换为 NaN（而非 epsilon 替换，避免人为偏移结果）
    mn_safe = mn.where(mn.abs() > EPSILON, np.nan)
    cn_safe = cn.where(cn.abs() > EPSILON, np.nan)

    ratio_mut = ms / mn_safe
    ratio_con = cs / cn_safe

    ratio_con_safe = ratio_con.where(ratio_con.abs() > EPSILON, np.nan)
    rf = ratio_mut / ratio_con_safe

    if feature in LOWER_IS_BETTER:
        rf_safe = rf.where(rf.abs() > EPSILON, np.nan)
        rf = 1.0 / rf_safe

    return rf


def main():
    parser = argparse.ArgumentParser(description='Calculate RF Matrix for Selected Features')
    parser.add_argument('--input_features', required=True,
                        help='Path to all_features_raw_50plus.csv')
    parser.add_argument('--input_selected', required=True,
                        help='Path to selected_features_list.csv (must have a "Feature" column)')
    parser.add_argument('--output_dir', required=True, help='Output directory')
    args = parser.parse_args()

    # ---- 加载数据 ----
    print(f"Loading features from {args.input_features}...")
    df = pd.read_csv(args.input_features)

    print(f"Loading selected features list from {args.input_selected}...")
    df_selected = pd.read_csv(args.input_selected)
    selected_features = df_selected['Feature'].tolist()

    if 'gene' not in df.columns:
        print("Error: 'gene' column not found in input features.")
        return

    gene_list = df['gene'].tolist()
    print(f"Genes: {len(gene_list)}, Features to process: {len(selected_features)}")

    # ---- 向量化计算 RF ----
    rf_data = {}
    skipped = []

    for feature in selected_features:
        ms_col = f'mut_stress_{feature}'
        mn_col = f'mut_nonstress_{feature}'
        cs_col = f'con_stress_{feature}'
        cn_col = f'con_nonstress_{feature}'

        missing = [c for c in [ms_col, mn_col, cs_col, cn_col] if c not in df.columns]
        if missing:
            print(f"  Warning: Feature '{feature}' missing columns {missing}. Skipping.")
            skipped.append(feature)
            continue

        rf_series = _safe_rf(
            df[ms_col], df[mn_col],
            df[cs_col], df[cn_col],
            feature
        )
        rf_data[feature] = rf_series.values

    if skipped:
        print(f"\nSkipped {len(skipped)} features due to missing columns: {skipped}")

    # ---- 构建结果 DataFrame ----
    df_rf = pd.DataFrame(rf_data, index=gene_list)

    # 删除所有特征均为 NaN 的基因行
    before = len(df_rf)
    df_rf = df_rf.dropna(how='all')
    after = len(df_rf)
    if before > after:
        print(f"Dropped {before - after} genes with all-NaN RF values.")

    # ---- 保存 ----
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    output_csv = output_dir / 'rf_matrix_final.csv'
    df_rf.to_csv(output_csv)

    # 简单统计
    nan_frac = df_rf.isna().mean()
    high_nan = nan_frac[nan_frac > 0.3]
    if not high_nan.empty:
        print("\nWarning: Features with >30% NaN RF values (consider removing):")
        for feat, frac in high_nan.items():
            print(f"  {feat}: {frac:.1%} NaN")

    print(f"\nRF Matrix saved to: {output_csv}")
    print(f"Matrix shape: {df_rf.shape}  (genes × features)")


if __name__ == "__main__":
    main()
