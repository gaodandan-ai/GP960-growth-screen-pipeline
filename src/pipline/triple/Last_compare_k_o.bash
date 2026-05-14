#!/bin/bash
# 自动定位路径
SCRIPT_PATH=$(readlink -f "$0")
PIPELINE_DIR=$(dirname "$SCRIPT_PATH")
BASE_DIR=$(dirname $(dirname $(dirname "$PIPELINE_DIR")))

# 定义结果目录
K_RESULT_DIR="${BASE_DIR}/data/results/triple_replicate/thermal_k/02.cleaned_data"
O_RESULT_DIR="${BASE_DIR}/data/results/triple_replicate/thermal_TF_o/02.cleaned_data"
OUTPUT_DIR="${BASE_DIR}/data/results/triple_replicate/comparison_k_o"
MAPPING_FILE="${BASE_DIR}/data/raw/triple_replicate/plate_gene_mapping.xlsx"

# 设置最大时间限制（小时）
MAX_TIME=25

# 运行对比绘图脚本 (分条件汇总)
python ${BASE_DIR}/src/data_processing/comparison_k_o/01.plot_combined_curves.py \
    --k_dir ${K_RESULT_DIR} \
    --o_dir ${O_RESULT_DIR} \
    --output_dir ${OUTPUT_DIR} \
    --time_max ${MAX_TIME} \
    --mapping_file ${MAPPING_FILE}

# 运行对比绘图脚本 (同基因跨条件对比)
python ${BASE_DIR}/src/data_processing/comparison_k_o/02.plot_condition_comparison.py \
    --k_dir ${K_RESULT_DIR} \
    --o_dir ${O_RESULT_DIR} \
    --output_dir ${OUTPUT_DIR} \
    --time_max ${MAX_TIME} \
    --mapping_file ${MAPPING_FILE}
