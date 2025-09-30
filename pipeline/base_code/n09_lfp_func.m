%%% Base function of Micros Pipeline, Step 9: LFP
%%% This function saves a copy of downsampled data
%%% for LFP to be used in different analyses.

function n09_lfp_func(varargin)
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
clean_directory = fullfile(data_directory, subject, folder, 'clean', sprintf('Bank%s', bank));
lfp_directory = strrep(clean_directory, 'clean', 'lfp');
if ~exist(lfp_directory, 'dir')
    mkdir(lfp_directory)
end

channel_files = dir(fullfile(clean_directory, 'NS6*.mat'));
channel_files = {channel_files.name};

for idx = 1:length(channel_files)
    file_name = channel_files{idx};
    clean_file = fullfile(clean_directory, file_name);
    
    load(clean_file, 'data')
    
    lfp = downsample(data, 15);  %%% Not rescaling to conserve resolution
    sampling_rate = 2000;
    
    file_name = strrep(file_name, 'NS6', 'NS3');
    save(fullfile(lfp_directory, file_name), 'lfp', 'sampling_rate');
end
end
