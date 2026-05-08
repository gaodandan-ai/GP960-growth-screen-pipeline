#!/bin/bash
# step1: pp_rawdata.py
python ../data_processing/triple_replicate/01.pp_rawdata.py \
    --rawdata_dir ../../data/raw/triple_replicate/thermal_k \
    --plate_file ../../data/raw/triple_replicate/knockout_plate_gene_mapping.xlsx \
    --result_dir ../../data/results/triple_replicate/thermal_k/01.ppraw_data \
    --c1_label nonstress \
    --c2_label stress

## step2: data_cleaning.py
python ../data_processing/triple_replicate/02.data_cleaning.py \
    --input_dir ../../data/results/triple_replicate/thermal_k/01.ppraw_data \
    --output_dir ../../data/results/triple_replicate/thermal_k/02.cleaned_data \
    --mapping_file ../../data/raw/triple_replicate/knockout_plate_gene_mapping.xlsx \
    --threshold 0

## step3: calculate_fitness.py
python ../data_processing/triple_replicate/03.calculate_fitness.py \
    --input_dir ../../data/results/triple_replicate/thermal_k/02.cleaned_data \
    --output_dir ../../data/results/triple_replicate/thermal_k/03.stress_tolerance_analysis \
    --rf_threshold_stress 0 \
    --comprehensive_rf_threshold_stress 0 \
    --rf_threshold_nonstress 0.8 \
    --comprehensive_rf_threshold_nonstress 0.8

## step4: plot_gene_od.py
python ../data_processing/triple_replicate/04.plot_gene_od.py \
    --input_dir ../../data/results/triple_replicate/thermal_k/02.cleaned_data \
    --output_dir ../../data/results/triple_replicate/thermal_k/04.gene_growth_curves \
    --analysis_dir ../../data/results/triple_replicate/thermal_k/03.stress_tolerance_analysis \
    --time_max 48.0

## step5: calculate_full_features.py
python ../data_processing/triple_replicate/05.calculate_full_features.py \
    --input_dir ../../data/results/triple_replicate/thermal_k/02.cleaned_data \
    --output_dir ../../data/results/triple_replicate/thermal_k/05.full_features

## step6: plot_raw_correlation_heatmap_exact.py
python ../data_processing/triple_replicate/06.plot_raw_correlation_heatmap_exact.py \
    --input ../../data/results/triple_replicate/thermal_k/05.full_features/all_features_raw_50plus.csv \
    --output_dir ../../data/results/triple_replicate/thermal_k/06.plot_feature_selection

## step7: calculate_rf_matrix.py
python ../data_processing/triple_replicate/07.calculate_rf_matrix.py \
    --input_features ../../data/results/triple_replicate/thermal_k/05.full_features/all_features_raw_50plus.csv \
    --input_selected ../../data/results/triple_replicate/thermal_k/06.plot_feature_selection/selected_features_list.csv \
    --output_dir ../../data/results/triple_replicate/thermal_k/07.select_features_RF

## step8: screen_top5_per_feature.py
python ../data_processing/triple_replicate/08.screen_top5_per_feature.py \
    --input ../../data/results/triple_replicate/thermal_k/07.select_features_RF/rf_matrix_final.csv \
    --output_dir ../../data/results/triple_replicate/thermal_k/08.top5_strains \
    --top_n 5

## step9: plot_top_strains_curves.py
python ../data_processing/triple_replicate/09.plot_top_strains_curves.py \
    --input_dir ../../data/results/triple_replicate/thermal_k/02.cleaned_data \
    --top5_file ../../data/results/triple_replicate/thermal_k/08.top5_strains/top5_strains_per_feature.csv \
    --output_dir ../../data/results/triple_replicate/thermal_k/09.plot_top5_curves_by_feature

## step10: plot_dotplot.py
python ../data_processing/triple_replicate/10.plot_dotplot.py \
    --top5_file ../../data/results/triple_replicate/thermal_k/08.top5_strains/top5_strains_per_feature.csv \
    --output_dir ../../data/results/triple_replicate/thermal_k/10.dotplot
