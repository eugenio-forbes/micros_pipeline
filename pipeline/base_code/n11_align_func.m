%%% Base function of Micros Pipeline, Step 10: Behavioral Data Alignment
%%% This function will align behavioral data to the sync pulses delivered
%%% to EEG recording so that spike data analyses related to behavior can
%%% be performed.

function [error_flag,error_message] =  n11_align_func(varargin)
if isempty(varargin)                                               %%% To run manually edit values below
    root_directory = '/path/to/micros_pipeline/parent_directory';  %%% Root directory with pipeline and database folders
    subject = 'SC000';                                             %%% Subject code in SC000 format
    folder = 'yyyy-mm-dd_task-code_part1';                         %%% Folder in yyyy-mm-dd_task-code format, or yyyy-mm-dd_task-code_part1 format if more than 1 part. 
    sync_num = 129;                                                %%% Number of DC channel containing sync pulse recording
else                                                               %%% Otherwise this is the order they should be entered into function, following above format
    root_directory = varargin{1};
    subject = varargin{2};
    folder = varargin{3};
    sync_num = varargin{4};
end

data_directory = fullfile(root_directory, 'micros_database', subject, folder);
behavioral_directory = fullfile(data_directory, 'behavioral');
split_directory = fullfile(data_directory, 'split');
lfp_directory = fullfile(data_directory, 'lfp', 'Bank*');
combinato_directory = fullfile(data_directory, 'combinato_files');
report_directory = fullfile(data_directory, 'alignment');

if ~exist(report_directory, 'dir')
    mkdir(report_directory);
end

error_flag = 0;
error_message = '';

events_file = fullfile(behavioral_directory, 'events.mat');
if isfile(fullfile(behavioral_directory, 'eeg.eeglog.up'))
    sync_file = fullfile(behavioral_directory, 'eeg.eeglog.up');
elseif isfile(fullfile(behavioral_directory, 'pulses.txt'))
    sync_file = fullfile(behavioral_directory, 'pulses.txt');
else
    error_flag = 1;
    error_message = 'No behavioral sync pulses in folder.';
    return
end

[downsampled_rate, n_samples] = get_LFP_info(lfp_directory);
lfp_directory = fullfile(data_directory, 'lfp');

%%% Pulse matching
blackrock_pulses = get_blackrock_pulses(split_directory, sync_num); %%% Pulses are found in NS6_129.mat, for example
file_id = fopen(sync_file, 'r');
behavioral_pulses = cell2mat(textscan(file_id, '%f%*[^\n]', 'Delimiter', '\t')); fclose(file_id);

%%% Try different window sizes for good alignment (some cases fail)
successful_alignment = false;
attempt_number = 1;
window_sizes = [30, 50, 25, 75, 20, 100, 15, 200, 10, 300, 5, 500, 3];
while ~successful_alignment && attempt_number <= length(window_sizes)
    window_size = window_sizes(attempt_number);
    [behavioral_matched, blackrock_matched] = align_blackrock_pulses(behavioral_pulses, blackrock_pulses', window_size);
    blackrock_sampling_rate = round(((blackrock_matched(end) - blackrock_matched(1)) / (behavioral_matched(end) - behavioral_matched(1))) * 1000); %%%29998; %%%MT added: calc sample rate from pulses and tot time
    [aligned_events, report_stats, successful_alignment] = align_events(events_file, lfp_directory, behavioral_matched, blackrock_matched, blackrock_sampling_rate, downsampled_rate, n_samples, combinato_directory);
    attempt_number = attempt_number + 1;
end

if successful_alignment
    events = aligned_events;
    
    movefile(fullfile(behavioral_directory, 'events.mat'), fullfile(behavioral_directory, 'unaligned_events.mat'));
    save(fullfile(behavioral_directory, 'events.mat'), 'events');
    
    report_file = fullfile(report_directory, 'report_stats.mat');
    save(report_file, 'report_stats')
    
    report_file = strrep(report_file, '.mat', '.txt')
    
    file_id = fopen(report_file, 'w');
    fprintf(file_id, 'Max. deviation  : %.02f ms\n', report_stats.max_deviation);
    fprintf(file_id, 'Med. deviation  : %.02f ms\n', report_stats.median_deviation);
    fprintf(file_id, '95th percentile : %.02f ms\n', report_stats.t95th_pctile);
    fprintf(file_id, '99th percentile : %.02f ms\n', report_stats.t99th_pctile);
    fprintf(file_id, 'R-squared       : %.05f \n', report_stats.R_squared);
    fprintf(file_id, 'Slope           : %.05f \n', report_stats.tSlope);
    fprintf(file_id, 'Pulse_range     : %.02f minutes\n', report_stats.pulse_range);
    fclose(file_id);
    
else
    error_flag = 1;
    error_message = 'Unsuccessful alignment of events with pulses and several window sizes.';
    report_file = fullfile(report_directory, 'unsuccessful_alignment');
    file_id = fopen(report_file, 'w');
    fclose(file_id);
end

end


function [downsampled_rate, n_samples] = get_LFP_info(lfp_directory)
channels = dir(fullfile(lfp_directory, 'NS3_*.mat'));
channels = fullfile({channels.folder}, {channels.name});
channel = channels{1};
load(channel, 'lfp', 'sampling_rate')
downsampled_rate = sampling_rate;
n_samples = length(lfp);
end


function blackrock_pulses = get_blackrock_pulses(split_directory, sync_num)
load(fullfile(split_directory, sprintf('NS6_%03d.mat', sync_num)), 'data');
data(data<0) = 0;
normalizedpulse = data/max(data);
dPulse = [0 diff(normalizedpulse)];
ddPulseShift = [diff(dPulse) 0];
dPulse(ddPulseShift==-2) = 0;
triggerThresh = 0.5;
%%% detect rising edge
blackrock_pulses = find(dPulse>triggerThresh);
end


function [aligned_events, report_stats, successful_alignment] = align_events(events_file, lfp_directory, behavioral_matched, blackrock_matched, blackrock_sampling_rate, downsampled_rate, n_samples, times_directory) %%%br_ms actually in samples, not ms

%%% allocate
b = zeros(2, 1);

%%% load events
aligned_events = load(events_file);
aligned_events = aligned_events.events;

%%% add blank fields
[aligned_events.lfpfile] = deal('');
[aligned_events.lfpoffset] = deal(NaN);
[aligned_events.timesoffset] = deal(NaN);

%%% behavioral times
event_times = [aligned_events.mstime]';

bfix = behavioral_matched(1);
[b(:), ~, ~, ~, stats] = regress(blackrock_matched, [ones(length(behavioral_matched), 1) behavioral_matched-bfix]);
b(1) = b(1) - bfix*b(2);
    
%%% calc max deviation
act=[ones(length(behavioral_matched), 1) behavioral_matched]*b(:);
maxdev = max(abs(act - blackrock_matched));
meddev = median(abs(act - blackrock_matched));
rawDevs = act - blackrock_matched;
    
%%% report stats
report_stats.max_deviation = maxdev*1000/blackrock_sampling_rate;
report_stats.median_deviation = meddev*1000/blackrock_sampling_rate;
report_stats.t95th_pctile = prctile(rawDevs, 95)*1000/blackrock_sampling_rate;
report_stats.t99th_pctile = prctile(rawDevs, 99)*1000/blackrock_sampling_rate;
report_stats.r_squared = stats(1);
report_stats.slope = b(2);
report_stats.pulse_range = range(behavioral_matched)/(1000*60);

successful_alignment = report_stats.max_deviation<2; %sometimes it's a little bit over 1 ms
   
%%% get start and stop in ms of the lfp file
eeg_start_ms = round((1 - b(1))/b(2));
eeg_stop_ms  = round((n_samples*(blackrock_sampling_rate/downsampled_rate) - b(1)/b(2)));
    
%%% convert blackrock sample to ms (for spikes) and to samplenum for
%%% lfpoffset. brEvents_times is really the sample number in the original
%%% (not downsampled) black rock file
brEvents_times = double([ones(length(event_times), 1) event_times])*b(:);
    
%%% convert to ms since start of recording (timesoffset)
brEvents_ms = brEvents_times / (blackrock_sampling_rate/1e3);
    
%%% convert to samples at 2000 hz (lfpoffset)
brEvents_samps = round(brEvents_times / (blackrock_sampling_rate/downsampled_rate));
    
%%% loop over each events
for e = 1:length(aligned_events)
    %%% check if this event is with the calculated bounds of the current
    %%% lfp file. Only add the timing info if it is.
    if aligned_events(e).mstime >= eeg_start_ms && aligned_events(e).mstime <= eeg_stop_ms
        aligned_events(e).NS3_lfpfile = lfp_directory;
        %%% in downsampled lfp samples
        aligned_events(e).NS3_lfpoffset_samples = brEvents_samps(e);
        aligned_events(e).timesfile = times_directory;
        %%% in ms
        aligned_events(e).timesoffset_ms = brEvents_ms(e);
    end 
end
end
