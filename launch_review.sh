#!/bin/bash

ROOT_DIR="/path/to/micros_pipeline/parent_directory" #Change to folder containing pipeline and database
TYPE=$1 #In command terminal must enter type: 'manual', 'task', 'subject'
FOLDER=$2 #In command terminal must enter folder: 'all' or task name for 'manual' and 'task' types, subject code for 'subject' type

# Run MATLAB code without display to change combinato jobs list.
module load matlab
matlab -nodesktop -nodisplay -singleCompThread -r "launch_review('$ROOT_DIR', '$TYPE', '$FOLDER'); exit;"

# Load anaconda and source conda.sh file so that conda commands can be executed here
module load python/3.6.4-anaconda
source /cm/shared/apps/python/3.6.4-anaconda/etc/profile.d/conda.sh
export PATH="${PATH}:/path/to/micros_pipeline/pipeline/base_code/combinato"
export PYTHONPATH="${PYTHONPATH}:/path/to/micros_pipeline/pipeline/base_code/combinato"

# Check if combinato conda environment exists from list and activate it
if $(conda env list | grep -q "combinato")
then
    conda activate combinato
else                          #If not then install dependencies from .yml and add paths to .bashrc

    conda env create -f /path/to/micros_pipeline/pipeline/base_code/combinato.yml
    echo 'export PATH="${PATH}:/path/to/micros_pipeline/pipeline/base_code/combinato"' >> ~/.bashrc
    echo 'export PYTHONPATH="${PYTHONPATH}:/path/to/micros_pipeline/pipeline/base_code/combinato"' >> ~/.bashrc
    echo 'source /cm/shared/apps/python/3.6.4-anaconda/etc/profile.d/conda.sh' >> ~/.bashrc
    conda activate combinato
fi

# Change directory to folder containing combinato_files and jobs lists and open combinato gui there
cd $ROOT_DIR
cd micros_pipeline/process_files/n07_manual/combinato_files
css-gui