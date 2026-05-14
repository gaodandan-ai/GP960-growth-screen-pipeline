# Set base directory relative to the script location
SCRIPT_PATH=$(readlink -f "$0")
PIPELINE_DIR=$(dirname "$SCRIPT_PATH")
BASE_DIR=$(dirname $(dirname $(dirname "$PIPELINE_DIR")))

# Set directories relative to the repository root
RAW_DIR="${BASE_DIR}/data/raw/single_replicate/thermal_overpress_all"
MAPPING_FILE="${RAW_DIR}/plate_gene_mapping.xlsx"
RESULT_BASE="${BASE_DIR}/data/results/single_replicate/thermal_overpress_all"
SCRIPT_DIR="${BASE_DIR}/src/data_processing/single_replicate_thermal_all"

# Ensure result directory exists
mkdir -p ${RESULT_BASE}

# Settings
MAX_TIME=24.0 

# step1: pp_rawdata.py
python ${SCRIPT_DIR}/01.pp_rawdata.py \
     --rawdata_dir ${RAW_DIR} \
     --plate_file ${MAPPING_FILE} \
     --result_dir ${RESULT_BASE}/01.ppraw_data \
     --c1_label stress \
     --c2_label nonstress

# step2: data_cleaning.py
python ${SCRIPT_DIR}/02.data_cleaning.py \
     --input_dir ${RESULT_BASE}/01.ppraw_data \
     --output_dir ${RESULT_BASE}/02.cleaned_data \
     --mapping_file ${MAPPING_FILE} \
     --threshold 1

# step3: calculate_fitness.py
python ${SCRIPT_DIR}/03.calculate_fitness.py \
    --input_dir ${RESULT_BASE}/02.cleaned_data \
    --output_dir ${RESULT_BASE}/03.stress_tolerance_analysis \
    --rf_threshold_stress 1.05 \
    --comprehensive_rf_threshold_stress 1.05 \
    --rf_threshold_nonstress 0.95 \
    --comprehensive_rf_threshold_nonstress 0.95 \
    --time_max ${MAX_TIME}

# step4: plot_gene_od.py
python ${SCRIPT_DIR}/04.plot_gene_od.py \
    --input_dir ${RESULT_BASE}/02.cleaned_data \
    --output_dir ${RESULT_BASE}/04.gene_growth_curves \
    --analysis_dir ${RESULT_BASE}/03.stress_tolerance_analysis \
    --time_max ${MAX_TIME} \
    --mapping_file ${MAPPING_FILE}

# step5: calculate_full_features.py
python ${SCRIPT_DIR}/05.calculate_full_features.py \
    --input_dir ${RESULT_BASE}/02.cleaned_data \
    --output_dir ${RESULT_BASE}/05.full_features \
    --time_max ${MAX_TIME}

# step6: plot_raw_correlation_heatmap_exact.py
python ${SCRIPT_DIR}/06.plot_raw_correlation_heatmap_exact.py \
    --input ${RESULT_BASE}/05.full_features/all_features_raw_50plus.csv \
    --output_dir ${RESULT_BASE}/06.plot_feature_selection

# step6: plot_raw_feature_scatter.py
python ${SCRIPT_DIR}/06.plot_raw_feature_scatter.py \
    --input ${RESULT_BASE}/05.full_features/all_features_raw_50plus.csv \
    --output_dir ${RESULT_BASE}/06.plot_feature_selection \
    --feat1 max_OD \
    --feat2 AUC

# step7: calculate_rf_matrix.py
python ${SCRIPT_DIR}/07.calculate_rf_matrix.py \
    --input_features ${RESULT_BASE}/05.full_features/all_features_raw_50plus.csv \
    --input_selected ${RESULT_BASE}/06.plot_feature_selection/selected_features_list.csv \
    --output_dir ${RESULT_BASE}/07.select_features_RF

# step8: screen_top5_per_feature.py
python ${SCRIPT_DIR}/08.screen_top5_per_feature.py \
    --input ${RESULT_BASE}/07.select_features_RF/rf_matrix_final.csv \
    --output_dir ${RESULT_BASE}/08.top5_strains \
    --top_n 5

# step9: plot_top_strains_curves.py
python ${SCRIPT_DIR}/09.plot_top_strains_curves.py \
    --input_dir ${RESULT_BASE}/02.cleaned_data \
    --top5_file ${RESULT_BASE}/08.top5_strains/top5_strains_per_feature.csv \
    --output_dir ${RESULT_BASE}/09.plot_top5_curves_by_feature \
    --time_max ${MAX_TIME} \
    --mapping_file ${MAPPING_FILE}

# step10: plot_dotplot.py
python ${SCRIPT_DIR}/10.plot_dotplot.py \
    --top5_file ${RESULT_BASE}/08.top5_strains/top5_strains_per_feature.csv \
    --output_dir ${RESULT_BASE}/10.dotplot
