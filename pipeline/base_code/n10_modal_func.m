%%% Base function of Micros Pipeline, Step 10: MODAL
%%% This function will perform oscillation detection
%%% on the LFP of every channel in a microelectrode,
%%% using MODAL function.

function n10_modal_func(varargin)
if isempty(varargin)                                               %%% To run manually edit values below
    root_directory = '/path/to/micros_pipeline/parent_directory';  %%% Root directory with pipeline and database folders
    subject = 'SC000';                                             %%% Subject code in SC000 format
    folder = 'yyyy-mm-dd_task-code_part1';                         %%% Folder in yyyy-mm-dd_task-code format, or yyyy-mm-dd_task-code_part1 format if more than 1 part. 
    bank = 'A';                                                    %%% Recording hardware bank character ('A', 'B', 'C', 'D')
else                                                               %%% Otherwise this is the order they should be entered into function, following above format
    root_directory = varargin{1};
    subject = varargin{2};
    folder = varargin{3};
    bank = varargin{4};
end
data_directory = fullfile(root_directory, 'micros_database');
lfp_directory = fullfile(data_directory, subject, folder, 'clean', sprintf('Bank%s', bank));
modal_directory = strrep(lfp_directory, 'lfp', 'modal');
if ~exist(modal_directory, 'dir')
    mkdir(modal_directory)
end

channel_files = dir(fullfile(lfp_directory, 'NS3*.mat'));
channel_files = {channel_files.name};

%%% Parameters for MODAL
frequencies = (2.^((8:40)/8));             %%% Frequencies from 2 to 32 Hz to find low frequency oscillations
modal_parameters.srate = 2000;             %%% Set this to the resampled rate of the LFPs: 2000 Hz
modal_parameters.wavefreqs = frequencies;
modal_parameters.crop_fs = 1;

for idx = 1:length(channel_files)
    file_name = channel_files{idx};
    lfp_file = fullfile(lfp_directory, file_name);
    load(lfp_file, 'lfp')
   
    [frequency_sliding, frequency_bands, band_power, band_phase] = MODAL(double(lfp), modal_parameters);
    
    peak_frequency = zeros(size(frequency_sliding, 1), 1);
    for jdx = 1:size(frequency_sliding, 1)
        peak_frequency(jdx, 1) = mode(round(frequency_sliding(jdx, :), 2));
    end
    
    save(fullfile(modal_directory, file_name), 'frequency_sliding', 'frequency_bands', 'band_power', 'band_phase', 'modal_parameters', 'peak_frequency');
end
end
