%%% Base function of Micros Pipeline, Step 4: Rescaling
%%% This function saves a copy of rescaled signal for
%%% every copy of a microelectrode. Rescaled by a factor
%%% of 4 so that the data is in uV. Additionally
%%% calculates standard deviation of signal.

function n04_rescale_func(varargin)
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

data_directory = fullfile(root_directory, 'micros_database', subject, folder, 'clean', sprintf('Bank%s', bank));
split_directory = fullfile(root_directory, 'micros_database', subject, folder, 'split');
save_directory = strrep(data_directory, 'clean', 'rescaled');
if ~exist(save_directory, 'dir')
    mkdir(save_directory);
end

file_info = load(fullfile(split_directory, 'file_length.mat'));
sampling_rate = file_info.sampling_rate;
file_length = file_info.file_length;

file_paths = dir(fullfile(data_directory, '*.mat'));
file_names = {file_paths.name};
file_names = file_names(~contains(file_names, 'SNR'));

filtered_data = zeros(length(file_names), file_length);
for idx = 1:length(file_names)
    this_file = file_names{idx};
    load(fullfile(data_directory, this_file), 'data')
    
    data = double(data);
    filtered_data(idx, :) = round(0.25 * data);
    data = int16(filtered_data(idx, :));
    
    save(fullfile(save_directory, this_file), 'data', 'sampling_rate');
    clear data
end

min_pass = 300; %Hz
max_pass = 3000; %Hz
order = 3;
[butterworth_b, butterworth_a] = butter(order, [min_pass, max_pass]/(sampling_rate/2), 'bandpass');

filtered_data = filtfilt(butterworth_b, butterworth_a, filtered_data')';    

noise_standard_deviations = median(abs(filtered_data), 2) ./ 0.6745; 

clear filtered_data

if ~iscolumn(file_names)
    file_names = file_names';
end

noise_info = table;
noise_info.file_name = file_names;
noise_info.noise_standard_deviation = noise_standard_deviations;

save(fullfile(save_directory, 'noise_info.mat'), 'noise_info');

end