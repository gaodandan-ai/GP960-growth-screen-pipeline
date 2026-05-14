#!/bin/bash
# 自动定位路径
SCRIPT_PATH=$(readlink -f "$0")
PIPELINE_DIR=$(dirname "$SCRIPT_PATH")
BASE_DIR=$(dirname $(dirname $(dirname "$PIPELINE_DIR")))
PROC_DIR="${BASE_DIR}/src/data_processing/triple_replicate"
RESULT_DIR="${BASE_DIR}/data/results/triple_replicate/thermal_TF_o"
RAW_DIR="${BASE_DIR}/data/raw/triple_replicate"

# 设置最大时间限制（小时）
MAX_TIME=30.0

# step1: pp_rawdata.py
python ${PROC_DIR}/01.pp_rawdata.py \
    --rawdata_dir ${RAW_DIR}/thermal_TF_o \
    --plate_file ${RAW_DIR}/TFoverpress_plate_gene_mapping.xlsx \
    --result_dir ${RESULT_DIR}/01.ppraw_data \
    --c1_label nonstress \
    --c2_label stress

## step2: data_cleaning.py
python ${PROC_DIR}/02.data_cleaning.py \
    --input_dir ${RESULT_DIR}/01.ppraw_data \
    --output_dir ${RESULT_DIR}/02.cleaned_data \
    --mapping_file ${RAW_DIR}/TFoverpress_plate_gene_mapping.xlsx \
    --threshold 1

## step3: calculate_fitness.py
python ${PROC_DIR}/03.calculate_fitness.py \
    --input_dir ${RESULT_DIR}/02.cleaned_data \
    --output_dir ${RESULT_DIR}/03.stress_tolerance_analysis \
    --rf_threshold_stress 0 \
    --comprehensive_rf_threshold_stress 0 \
    --rf_threshold_nonstress 0.8 \
    --comprehensive_rf_threshold_nonstress 0.8 \
    --time_max ${MAX_TIME}

## step4: plot_gene_od.py
python ${PROC_DIR}/04.plot_gene_od.py \
    --input_dir ${RESULT_DIR}/02.cleaned_data \
    --output_dir ${RESULT_DIR}/04.gene_growth_curves \
    --analysis_dir ${RESULT_DIR}/03.stress_tolerance_analysis \
    --time_max ${MAX_TIME} \
    --mapping_file ${RAW_DIR}/TFoverpress_plate_gene_mapping.xlsx

## step5: calculate_full_features.py
python ${PROC_DIR}/05.calculate_full_features.py \
    --input_dir ${RESULT_DIR}/02.cleaned_data \
    --output_dir ${RESULT_DIR}/05.full_features \
    --time_max ${MAX_TIME}

## step6: plot_raw_correlation_heatmap_exact.py
python ${PROC_DIR}/06.plot_raw_correlation_heatmap_exact.py \
    --input ${RESULT_DIR}/05.full_features/all_features_raw_50plus.csv \
    --output_dir ${RESULT_DIR}/06.plot_feature_selection

## step7: calculate_rf_matrix.py
python ${PROC_DIR}/07.calculate_rf_matrix.py \
    --input_features ${RESULT_DIR}/05.full_features/all_features_raw_50plus.csv \
    --input_selected ${RESULT_DIR}/06.plot_feature_selection/selected_features_list.csv \
    --output_dir ${RESULT_DIR}/07.select_features_RF

## step8: screen_top5_per_feature.py
python ${PROC_DIR}/08.screen_top5_per_feature.py \
    --input ${RESULT_DIR}/07.select_features_RF/rf_matrix_final.csv \
    --output_dir ${RESULT_DIR}/08.top5_strains \
    --top_n 5

## step9: plot_top_strains_curves.py
python ${PROC_DIR}/09.plot_top_strains_curves.py \
    --input_dir ${RESULT_DIR}/02.cleaned_data \
    --top5_file ${RESULT_DIR}/08.top5_strains/top5_strains_per_feature.csv \
    --output_dir ${RESULT_DIR}/09.plot_top5_curves_by_feature \
    --time_max ${MAX_TIME} \
    --mapping_file ${RAW_DIR}/TFoverpress_plate_gene_mapping.xlsx

## step10: plot_dotplot.py
python ${PROC_DIR}/10.plot_dotplot.py \
    --top5_file ${RESULT_DIR}/08.top5_strains/top5_strains_per_feature.csv \
    --output_dir ${RESULT_DIR}/10.dotplot
