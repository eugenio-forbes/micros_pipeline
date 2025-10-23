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
- Save spike data in lighter format.

## Example Plots
The key to identifying neural activity in microelectrode recordings is, first of all, ensuring that the signal is of good quality. If line noise power exceeds signal power, then most of the spikes that are extracted will be the peaks and troughs of the fundamental frequency and harmonics. Therefore, this project places significant emphasis in methods that would effectively eliminate noise, and provides visualization methods to assess this.

First of all, one may be able to reduce the noise that is present in all channels with a referencing method. Methods available in this project include referencing to a calculated common average reference, selecting the channel that had the lowest signal (spike bandwidth) power, selecting the channel that had the greatest line noise power, or not applying a software reference at all.

The figure below shows that selecting the channel with the lowest signal power as a reference was effective in reducing line noise power by close to 35dB.
<p align="center">
  <img src="/process_files/n02_reref/plots/2019-10-05_AR_BankB_signal_power.png" alt="Reref1" width="1080"/>
</p>

One has to be careful with the referencing method selected, since this may duplicate signal to other channels as seen in the example below. It could also increase line noise power in other channels when selecting the channel with the greatest amount of line noise.
<p align="center">
  <img src="/process_files/n02_reref/plots/2019-10-05_AR_BankA_signal_power.png" alt="Reref2" width="1080"/>
</p>

Another method for dealing with line noise would be filtering. Although combinato offers band pass filters of the signal, it does not handle harmonics. Zapline by [NoiseTools](http://audition.ens.fr/adc/NoiseTools/) may be effective at further reducing the power of line noise and harmonics, as shown below.
<p align="center">
  <img src="/process_files/n03_clean/plots/2019-10-05_AR_BankA.png" alt="Clean1" width="1080"/>
</p>

It may not always be possible to remove noise from very noisy channels.
<p align="center">
  <img src="/process_files/n03_clean/plots/2020-01-12_AR_BankB.png" alt="Clean2" width="1080"/>
</p>

Next, considering how timely it would be to review results to verify whether hundreds of thousands of spike clusters are actual single unit neural activity (SU) vs. noise/artifact (N), this project attempts to further automatize this process by comparing features and clustering quality metrics of labeled clusters. Below is an example of this comparison.
<p align="center">
  <img src="/quality_metrics/plots/neurons/2023-09-12_all_neurons_sua_vs_noise.png" alt="Quality1" width="1080"/>
</p>

Finally, one could also compare quality metrics between different signal acquisition devices, such as digital (D) vs. analog (A) headstages, to assess which device yields better results.
<p align="center">
  <img src="/quality_metrics/plots/neurons/2023-09-12_all_neurons_digital_vs_analog.png" alt="Quality2" width="1080"/>
</p>
