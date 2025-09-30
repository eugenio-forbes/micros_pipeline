# micros_pipeline
**This pipeline processes microelectrode recordings from Blackrock Neurotech to detect and cluster single unit neural activity.**

## Overview:
This project provides a Matlab-based pipeline that facilitates the processing of Blackrock Neurotech microelectrode recordings for use with spike clustering method [combinato](https://github.com/jniediek), with the objective of improving the yield and quality of single unit neural activity clusters. It is tailored for cognitive electrophysiologists working on microelectrode data research. The data produced can readily be used for cognitive electrophysiology analyses using other Matlab-based tools developed by Texas Computational Memory Lab.


## Features:
- **Multi-node Processing**: Pipeline can be executed from multiple nodes with 256GB of RAM to speed up processing of hundreds of microelectrode recordings.
- **Parallel Processing**: Memory usage from each pipeline step has been measured to determine number of workers that would maximize usage of RAM without exceeding it. 
- **Modular Structure**: New methods for referencing, denoising, and clustering can easily be integrated into the pipeline.
- **Spike Cluster Quality Metrics**: Measurement of quality metrics allows for easy comparison of the effectiveness of different clustering parameters. It also allows for more effective grading of spike clusters as single unit activity, multiuinit activity, or artifact; which would help in reducing the amount of time clustering results are manually reviewed.
- **Automatized Processing**: By integrating different code bases all processing can now be carried out by executing a single command in a terminal.


## Tech Stack:
- **Platform**: Linux.
- **Requirements**: Anaconda.
- **Languages**: Matlab.
- **Key Package**: Must install [combinato](https://github.com/jniediek).


## Installation:
1) Download repository into desired parent directory.
2) Install [combinato](https://github.com/jniediek) in micros_pipeline/pipeline/base_code folder. Anaconda environment for combinato is found in the same folder.
3) Edit paths for parent directory in all files.


## Usage:
1) Individual Blackrock Neurotech raw microelectrode recordings (.ns6) should be stored in folders with the following format: /parent_directory/micros_database/subject_code/yyyy-mm-dd_task/raw/microelectrode_recording.ns6

2) Execute micros_pipeline/update_micros.sh in terminal.


## Optional:
Edit parameters of functions in update_micros.sh

## Future Improvements:
- Add cleanup function to remove unnecessary intermediate files.
- Add more methods for referencing, denoising, and clustering.
- Plot single unit activity spikes back into raw signal to verify they are not noise.
- Accurate alignment would actually require correction of behavioral task computer time.
