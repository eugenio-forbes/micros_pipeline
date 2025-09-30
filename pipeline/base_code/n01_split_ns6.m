%%% Base function of Micros Pipeline, Step 1: Splitting Recording
%%% This function will take as inputs directory containing database folder,
%%% subject code, and task folder name, to localize a single .ns6 file
%%% containing the microelectrode recordings for corresponding experimental session.
%%% It will save the data of each channel in a microelectrode (8 channels each),
%%% in a unique file. It will also list the information of the recording
%%% in jacksheet and parameter files compatible with other lab functions.

function n01_split_ns6(varargin)
if isempty(varargin)                                               %%% To run manually edit values below
    root_directory = '/path/to/micros_pipeline/parent_directory';  %%% Root directory with pipeline and database folders
    subject = 'SC000';                                             %%% Subject code in SC000 format
    folder = 'yyyy-mm-dd_task-code_part1';                         %%% Folder in yyyy-mm-dd_task-code format, or yyyy-mm-dd_task-code_part1 format if more than 1 part.    
else                                                               %%% Otherwise this is the order they should be entered into function, following above format
    root_directory = varargin{1};
    subject = varargin{2};
    folder = varargin{3};
end

%%% List directories and add path to pipeline code
code_directory = fullfile(root_directory, 'micros_pipeline/base_code'); addpath(genpath(code_directory));
data_directory = fullfile(root_directory, 'micros_database', subject, folder);
raw_directory = fullfile(data_directory, 'raw');
split_directory = fullfile(data_directory, 'split');
if ~exist(split_directory, 'dir')
    mkdir(split_directory);
end

%%% List of file name parts that should not be included in splitting, the
%%% data format that .ns6 files are read and save in. The data should be
%%% divided by 4 to represent the true voltage. However dividing int16 data
%%% by 4 leads to a decrease in the resolution of the data, so for
%%% functions seeking true voltage, data should be loaded, converted to
%%% double and divided by 4 for use.

bad_tags = {'.txt', 'jacksheet', 'params'};
data_format = 'int16';
gain = 1;

%%% Find the .ns6 files (raw, 30kHz, unfiltered format in Blackrock).
%%% Make sure each folder only has one .ns6 file.
ns6_file = dir(fullfile(raw_directory, '*.ns6'));
if ~isempty(ns6_file)
ns6_file = fullfile({ns6_file.folder}, {ns6_file.name});
end
assert(length(ns6_file)==1, 'Expected 1 .ns6 file, found %d', length(ns6_file));

%%% Read the .ns6 file into a data structure 
EEG_recording = open_nsx(ns6_file{1});
if iscell(EEG_recording.Data) %%% this means data was fragmented into segments
    fragmented = 1;
else 
    fragmented = 0;  
end 

%%% Exclude non-channel files      
tag_labels = {EEG_recording.ElectrodesInfo.Label};
is_bad = ismember(tag_labels, bad_tags);
good_indices = find(~is_bad);

%%% Get info about recording date for creating file stems for jacksheet and
%%% params file. This could help verify that the right file was introduced
%%% into the right folder.
cell_date = num2cell(EEG_recording.MetaTags.DateTimeRaw([1:3, 5, 6]));
file_stem = fullfile(split_directory, sprintf('%s_%02d-%02d-%02d_%02d-%02d_ns6', subject, cell_date{:}));
jacksheet_file = strcat(file_stem, '.jacksheet.txt');
parameters_file = strcat(file_stem, '.params.txt');
length_file = fullfile(split_directory, 'file_length.mat');

sampling_rate = EEG_recording.MetaTags.SamplingFreq; %%% Sampling rate to be saved in .mat file for each channel.

%%% Changed the below code for saving split files and timestaps because
%%% rarely, even with only 1 NSP, the data will be divided in segments.
%%% Concatenating the data and adding time of all segments would solve the
%%% error. It goes into question whether there is a time interval with no
%%% recording between segments

%%% Save the recording data of each channel to an individual .mat file
for index = good_indices
    channel_file = strcat(sprintf('NS6_%03i', EEG_recording.ElectrodesInfo(index).ElectrodeID), '.mat');
    if fragmented == 0
        data = int16(EEG_recording.Data(index, :));
    elseif fragmented ==1
        data = [];
        for edx = 1:length(EEG_recording.Data)
            temp_data = EEG_recording.Data{edx};
            data = [channel_data, int16(temp_data(index, :))];
        end
    end
    save(fullfile(split_directory, channel_file), 'data', 'sampling_rate');
end
file_length = length(data);
save(length_file, 'file_length', 'sampling_rate');

%%% Write time stamps. The time stamps apply to all channels.
n_time_points = EEG_recording.MetaTags.DataPoints;
if length(n_time_points) > 1  %%% Sometimes there is more than one and this fixes error
    n_time_points = sum(n_time_points);
end
try 
    time_stamps = linspace(0, (n_time_points - 1) * 1e6 / sampling_rate, n_time_points);
catch
    time_stamps = linspace(single(0), single((n_time_points - 1) * 1e6 / sampling_rate), single(n_time_points));
end
time_stamps_file = fullfile(split_directory, 'NS6_time_stamps.mat');
save(time_stamps_file, 'time_stamps', 'n_time_points', 'sampling_rate');

%%% Write params.txt
print_parameters(EEG_recording, parameters_file, data_format, gain);

%%% Write jacksheet.txt
print_jacksheet(EEG_recording, good_indices, jacksheet_file);
end


%%% Function to print parameters to file
function print_parameters(EEG_recording, parameters_file, data_format, gain)
file_id = fopen(parameters_file, 'w', 'l');
fprintf(file_id, 'samplerate %.2f\ndataformat ''%s''\ngain %g\n', EEG_recording.MetaTags.SamplingFreq, data_format, gain);
fclose(file_id);
end


%%% Function to print jacksheet to file
function print_jacksheet(EEG_recording, good_indices, jacksheet_file)
file_id = fopen(jacksheet_file, 'w');
for index = good_indices
    channel_label = EEG_recording.ElectrodesInfo(index).Label;
    channel_number = EEG_recording.ElectrodesInfo(index).ElectrodeID;
    %%% Regular strtrim didn't work - weird characters at the end of the label
    fprintf(file_id, '%d %s\n', channel_number, regexp_trim(channel_label));
end
fclose(file_id);
end


%%% Function that uses regexp to trim non-alphanumeric characters
function new_string = regexp_trim(old_string)
new_string = regexprep(old_string, '[^a-zA-Z0-9]', '');
end