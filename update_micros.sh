#!/bin/bash

# Node start time, parent directory, and custom preferences for rereferencing,
# denoising, clustering

START_TIME=$(date +%Y%m%d%H%M%S)
ROOT_DIR="/path/to/micros_pipeline/parent_directory"   # Parent directory containing 'micros_pipeline' and 'micros_database' folders
REREF_METHOD="reref_func"                              # 
REREF_MODE="signal_power"
DENOISE_METHOD="zapline"
CLUSTER_METHOD="combinato_func"
CLUSTER_MODE="standard"

cd $ROOT_DIR
cd micros_pipeline
module load matlab

# Run MATLAB code without display
matlab -nodesktop -nodisplay -singleCompThread -r "update_micros('$START_TIME','$ROOT_DIR', '$REREF_METHOD', '$REREF_MODE', '$DENOISE_METHOD', '$CLUSTER_METHOD', '$CLUSTER_MODE'); exit;"
