function plot_quality_metrics(varargin)
if isempty(varargin)
    root_directory = '/path/to/micros_pipeline/parent_directory';
else
    root_directory = varargin{1};
end

%%% Declare directories and load progress table
data_directory      = fullfile(root_directory, 'micros_database');
quality_directory   = fullfile(root_directory, 'micros_pipeline/quality_metrics/plots');
progress_table_file = fullfile(data_directory, 'progress_table.mat');

load(progress_table_file, 'progress_table')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Spike Data Processing %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Filter table to collect all available spike data from the recordings that were successfully clustered and timed.
timed_recordings = progress_table(progress_table.spikes_timed, :);
[timed_recordings, channels, neurons, classes, unique_tasks] = process_spike_data(data_directory, timed_recordings);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Neurons Quality Metrics Analyses %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
%%% List of all variables to be statistically compared
neuron_variables = [{'n_spikes'; 'f_spikes'; 'n_groups'; 'n_classes'; 'stdNoise'}, ...
    {'p_spc'; 'p_1st'; 'p_2nd'; 'm_distances'; 'std_distances'}, ...
    {'mean_FR'; 'm_amp_uV'; 'std_amp_uV'; 'm_voltage'; 'std_voltage'}, ...
    {'mean_peak_uV'; 'peak_width'; 'p_sub3ms'; 'p_60Hz'; 'isi_SNR'}, ...
    {'dispersion'; 'm_mahal'; 'std_mahal'; 'isoD'; 'L_ratio'}, ...
    {'dunn_idx'; 'silhouette'; 'SNR_split'; 'SNR_reref'; 'SNR_clean'}];

%%% Corresponding plot titles
neuron_titles = [{'# of spikes in group'; 'fraction from all detections'; '# of groups in channel'; '# of classes in group'; 'STD of signal'}, ...
    {'proportion matched at SPC'; 'proportion matched at 1st template'; 'proportion matched at 2nd template'; 'mean distances'; 'STD distances'}, ...
    {'mean firing rate'; 'mean amplitude (uV)'; 'STD amplitude(uV)'; 'mean voltage(uV)'; 'STD voltage(uV)'}, ...
    {'mean waveform amplitude (uV)'; 'peak width (samples)'; 'proportion with ISI < 3ms'; 'proportion with ISIs of line noise'; 'ISI SNR'}, ...
    {'dispersion'; 'mean mahalanobis distance'; 'STD mahalanobis distance'; 'isolation distance'; 'L-ratio'}, ...
    {'Dunn index'; 'silhouette score'; 'SNR of raw signal'; 'SNR of rerefed signal'; 'SNR of denoised signal'}];

su_vs_noise_analysis(quality_directory, 'neurons', 'all-neurons', neuron_variables, neuron_titles, neurons);
su_vs_noise_analysis(quality_directory, 'neurons', 'digital', neuron_variables, neuron_titles, neurons(neurons.digital_headstage, :));
su_vs_noise_analysis(quality_directory, 'neurons', 'analog', neuron_variables, neuron_titles, neurons(~neurons.digital_headstage, :));
su_vs_noise_analysis(quality_directory, 'neurons', 'AR', neuron_variables, neuron_titles, neurons(neurons.task_ID == find(strcmp(unique_tasks, 'AR')), :));
su_vs_noise_analysis(quality_directory, 'neurons', 'SR', neuron_variables, neuron_titles, neurons(neurons.task_ID == find(strcmp(unique_tasks, 'SR')), :));

digital_vs_analog_analysis(quality_directory, 'neurons', 'all-neurons', neuron_variables, neuron_titles, neurons);
digital_vs_analog_analysis(quality_directory, 'neurons', 'single-units', neuron_variables, neuron_titles, neurons(neurons.is_SUA, :));
digital_vs_analog_analysis(quality_directory, 'neurons', 'noise-units', neuron_variables, neuron_titles, neurons(neurons.is_noise, :));
digital_vs_analog_analysis(quality_directory, 'neurons', 'AR', neuron_variables, neuron_titles, neurons(neurons.task_ID == find(strcmp(unique_tasks, 'AR')), :));
digital_vs_analog_analysis(quality_directory, 'neurons', 'SR', neuron_variables, neuron_titles, neurons(neurons.task_ID == find(strcmp(unique_tasks, 'SR')), :));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Classes Quality Metrics Analyses %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
%%% List of all variables to be statistically compared
class_variables = [{'n_spikes'; 'f_spikes'; 'n_groups'; 'n_classes'; 'stdNoise'}, ...
    {'p_spc'; 'p_1st'; 'p_2nd'; 'm_distances'; 'std_distances'}, ...
    {'L_ratio'; 'm_amp_uV'; 'std_amp_uV'; 'm_voltage'; 'std_voltage'}, ...
    {'m_dispersion'; 'm_mahal'; 'std_mahal'; 'intragroup_isoD'; 'extragroup_isoD'}, ...
    {'p_sub3ms'; 'p_60Hz'; 'SNR_split'; 'SNR_reref'; 'SNR_clean'}];
%%% Corresponding plot titles
class_titles = [{'# of spikes in class'; 'fraction from all detections'; '# of groups in channel'; '# of classes in group'; 'STD of signal'}, ...
    {'proportion matched at SPC'; 'proportion matched at 1st template'; 'proportion matched at 2nd template'; 'mean distances'; 'STD distances'}, ...
    {'L-ratio'; 'mean amplitude (uV)'; 'STD amplitude(uV)'; 'mean voltage(uV)'; 'STD voltage(uV)'}, ...
    {'dispersion'; 'mean mahalanobis distance'; 'std mahalanobis distance'; 'intragroup isolation distance'; 'extragroup isolation distance'}, ...
    {'proportion with ISI < 3ms'; 'proportion with ISIs of line noise'; 'SNR of raw signal'; 'SNR of rerefed signal'; 'SNR of denoised signal'}];

su_vs_noise_analysis(quality_directory, 'classes', 'all-classes', class_variables, class_titles, classes);
su_vs_noise_analysis(quality_directory, 'classes', 'digital', class_variables, class_titles, classes(classes.digital_headstage, :));
su_vs_noise_analysis(quality_directory, 'classes', 'analog', class_variables, class_titles, classes(~classes.digital_headstage, :));
su_vs_noise_analysis(quality_directory, 'classes', 'AR', class_variables, class_titles, classes(classes.task_ID == find(strcmp(unique_tasks, 'AR')), :));
su_vs_noise_analysis(quality_directory, 'classes', 'SR', class_variables, class_titles, classes(classes.task_ID == find(strcmp(unique_tasks, 'SR')), :));

digital_vs_analog_analysis(quality_directory, 'classes', 'all-classes', class_variables, class_titles, classes);
digital_vs_analog_analysis(quality_directory, 'classes', 'single-units', class_variables, class_titles, classes(classes.is_SUA, :));
digital_vs_analog_analysis(quality_directory, 'classes', 'noise-units', class_variables, class_titles, classes(classes.is_noise, :));
digital_vs_analog_analysis(quality_directory, 'classes', 'AR', class_variables, class_titles, classes(classes.task_ID == find(strcmp(unique_tasks, 'AR')), :));
digital_vs_analog_analysis(quality_directory, 'classes', 'SR', class_variables, class_titles, classes(classes.task_ID == find(strcmp(unique_tasks, 'SR')), :));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Channels Quality Metrics Analyses %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

channels.duration_seconds  = channels.file_length ./ channels.sr;
channels.detection_rate    = channels.n_detections ./ channels.duration_seconds;
channels.spikes_per_group  = channels.n_spikes ./ channels.n_groups;
channels.spikes_per_class  = channels.n_spikes ./ channels.n_classes;
channels.classes_per_group = channels.n_classes ./ channels.n_groups;
channels.p_spc             = channels.n_spc ./ channels.n_spikes;
channels.p_1st             = channels.n_1st ./ channels.n_spikes;
channels.p_2nd             = channels.n_2nd ./ channels.n_spikes;
channels.p_artifact        = channels.n_bad_spikes ./ channels.n_detections;

%%% List of all variables to be statistically compared
channel_ys = {'n_sua', 'n_noise'};
channel_ys_labels = {'# of SU', '# of MU', '# of N'};

channel_variables = [{'stdNoise'; 'detection_rate'; 'n_detections'; 'duration_seconds'; 'spikes_per_group'}, ...
    {'spikes_per_class'; 'classes_per_group'; 'p_spc'; 'p_1st'; 'p_2nd'}, ...
    {'p_artifact'; 'isi_SNR'; 'DBI'; 'CHI'; 'm_silhouette'}, ...
    {'SNR_split'; 'SNR_reref'; 'SNR_clean'; 'n_groups'; 'n_classes'}];

%%% Corresponding plot titles
channel_titles = [{'STD of signal'; 'detection rate'; '# of detections'; 'recording duration(s)'; '# spikes per group'}, ...
    {'#spikes per class'; '#classes per group'; 'proportion matched at SPC'; 'proportion matched at 1st template'; 'proportion matched at 2nd template'}, ...
    {'artifact spike fraction'; 'mean ISI SNR'; 'Davies-Bouldin Index'; 'Calinski-Harabasz Index'; 'mean silhouette score'}, ...
    {'SNR split'; 'SNR rerefed'; 'SNR denoised'; '# of groups'; '# of classes'}];

digital_vs_analog_analysis(quality_directory, 'channels', 'all-channels', channel_variables, channel_titles, channels);
for idx = 1:length(channel_ys)
    channel_y = channel_ys{idx};
    count_quality_correlation(quality_directory, 'channels', 'all-channels', channel_y, channel_variables, channel_titles, channels, 'electrodes');
    count_quality_correlation(quality_directory, 'channels', 'all-channels', channel_y, channel_variables, channel_titles, channels, 'subjects');
    count_quality_correlation(quality_directory, 'channels', 'all_channels', channel_y, channel_variables, channel_titles, channels, 'tasks');
end
digital_vs_analog_counts(quality_directory, 'channels', 'all-channels', channel_ys, channel_ys_labels, channels);

digital_vs_analog_analysis(quality_directory, 'channels', 'AR', channel_variables, channel_titles, channels(channels.task_ID == find(strcmp(unique_tasks, 'AR')), :));
for idx = 1:length(channel_ys)
    channel_y = channel_ys{idx};
    count_quality_correlation(quality_directory, 'channels', 'AR', channel_y, channel_variables, channel_titles, channels(channels.task_ID == find(strcmp(unique_tasks, 'AR')), :), 'electrodes');
    count_quality_correlation(quality_directory, 'channels', 'AR', channel_y, channel_variables, channel_titles, channels(channels.task_ID == find(strcmp(unique_tasks, 'AR')), :), 'subjects');
    count_quality_correlation(quality_directory, 'channels', 'AR', channel_y, channel_variables, channel_titles, channels(channels.task_ID == find(strcmp(unique_tasks, 'AR')), :), 'tasks');
end
digital_vs_analog_counts(quality_directory, 'channels', 'AR', channel_ys, channel_ys_labels, channels(channels.task_ID == find(strcmp(unique_tasks, 'AR')), :));

digital_vs_analog_analysis(quality_directory, 'channels', 'SR', channel_variables, channel_titles, channels(channels.task_ID == find(strcmp(unique_tasks, 'SR')), :));
for idx = 1:length(channel_ys)
    channel_y = channel_ys{idx};
    count_quality_correlation(quality_directory, 'channels', 'SR', channel_y, channel_variables, channel_titles, channels(channels.task_ID == find(strcmp(unique_tasks, 'SR')), :), 'electrodes');
    count_quality_correlation(quality_directory, 'channels', 'SR', channel_y, channel_variables, channel_titles, channels(channels.task_ID == find(strcmp(unique_tasks, 'SR')), :), 'subjects');
    count_quality_correlation(quality_directory, 'channels', 'SR', channel_y, channel_variables, channel_titles, channels(channels.task_ID == find(strcmp(unique_tasks, 'SR')), :), 'tasks');
end
digital_vs_analog_counts(quality_directory, 'channels', 'SR', channel_ys, channel_ys_labels, channels(channels.task_ID == find(strcmp(unique_tasks, 'SR')), :));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Electrodes Quality Metrics Analyses %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

electrode_variables = {'SNR_split', 'SNR_reref', 'SNR_clean', 'isi_SNR', 'su_counts', 'mu_counts', 'noise_counts'};
electrode_titles = {'SNR split', 'SNR reref', 'SNR clean', 'ISI SNR', '# of single units', '# of multi units', '# of noise units'};

electrode_quality_plots(quality_directory, 'electrodes', 'all-electrodes', electrode_variables, electrode_titles, timed_recordings, 'electrodes');
electrode_quality_plots(quality_directory, 'electrodes', 'all-electrodes', electrode_variables, electrode_titles, timed_recordings, 'subjects');
electrode_quality_plots(quality_directory, 'electrodes', 'all-electrodes', electrode_variables, electrode_titles, timed_recordings, 'digital');

electrode_quality_plots(quality_directory, 'electrodes', 'AR', electrode_variables, electrode_titles, timed_recordings(timed_recordings.task_ID == find(strcmp(unique_tasks, 'AR')), :), 'electrodes');
electrode_quality_plots(quality_directory, 'electrodes', 'AR', electrode_variables, electrode_titles, timed_recordings(timed_recordings.task_ID == find(strcmp(unique_tasks, 'AR')), :), 'subjects');
electrode_quality_plots(quality_directory, 'electrodes', 'AR', electrode_variables, electrode_titles, timed_recordings(timed_recordings.task_ID == find(strcmp(unique_tasks, 'AR')), :), 'digital');

electrode_quality_plots(quality_directory, 'electrodes', 'SR', electrode_variables, electrode_titles, timed_recordings(timed_recordings.task_ID == find(strcmp(unique_tasks, 'SR')), :), 'electrodes');
electrode_quality_plots(quality_directory, 'electrodes', 'SR', electrode_variables, electrode_titles, timed_recordings(timed_recordings.task_ID == find(strcmp(unique_tasks, 'SR')), :), 'subjects');
electrode_quality_plots(quality_directory, 'electrodes', 'SR', electrode_variables, electrode_titles, timed_recordings(timed_recordings.task_ID == find(strcmp(unique_tasks, 'SR')), :), 'digital');

end


function [timed_recordings, all_channels, all_neurons, all_classes, unique_tasks] = process_spike_data(data_directory, timed_recordings)

timed_recordings.folder = repelem({''}, height(timed_recordings), 1);
unique_channels = [];
bad_recordings = false(height(timed_recordings), 1);

for idx = 1:height(timed_recordings)
    
    %%% Get all folder names from date, task, and part
    subject = timed_recordings(idx, :).subject{:};
    date    = timed_recordings(idx, :).date{:};
    task    = timed_recordings(idx, :).task{:};
    part    = timed_recordings(idx, :).part;
    
    folder = strcat(date, '_', task);
    
    if part > 0
        folder = strcat(folder, '_', sprintf('part%d', part));
    end
    
    timed_recordings(idx, :).folder = {folder};
    
    %%% Get min and max channel numbers in bank
    channel_numbers = timed_recordings(idx, :).channels{:};
    
    %%% Avoid cases where jacksheet labels did not allow for accurate identification of a bank's channel numbers
    %%% These would have to be corrected manually
    if ~isempty(channel_numbers) && sum(channel_numbers) > 0
        
        %%% Get all channel numbers
        channel_numbers = channel_numbers(1):channel_numbers(2);
        channel_numbers = arrayfun(@(x) num2str(x), channel_numbers', 'UniformOutput', false);
        
        %%% Concatenate string to subject code to be able to identify all unique channels
        these_channels = strcat(subject, '_', channel_numbers);
        unique_channels = [unique_channels;these_channels];
        
    else
        bad_recordings(idx) = true;
    end
    
end

%%% Filter out recordings with innacurate channel number identification
timed_recordings(bad_recordings, :) = [];

%%% Get unique subjects, sessions, electrodes, and channels IDd
[~, ~, subject_IDs]         = unique(timed_recordings.subject, 'stable');
[~, ~, electrode_IDs]       = unique(strcat(timed_recordings.subject, '_', timed_recordings.bank), 'stable');
[~, ~, session_IDs]         = unique(strcat(timed_recordings.subject, '_', timed_recordings.folder), 'stable');
[unique_tasks, ~, task_IDs] = unique(timed_recordings.task, 'stable');
unique_channels             = unique(unique_channels, 'stable');

timed_recordings.subject_ID   = subject_IDs;
timed_recordings.electrode_ID = electrode_IDs;
timed_recordings.session_ID   = session_IDs;
timed_recordings.task_ID      = task_IDs;
    
%%% Declare empty tables to collect all channels, neurons, and classes
all_channels = table;
all_neurons = table;
all_classes = table;

%%% Loop through recordings table to load all channels, neurons, and classes.
%%% Associate them to their respective unique IDs

for idx = 1:height(timed_recordings)

    this_recording = timed_recordings(idx, :);
    
    %%% Two sets of channel numbers:
    
    %%% One for indexing SNR of split and rerefed data (rerefed contains one NaN)
    this_recording_channels = this_recording.channels{:};
    this_recording_channels = this_recording_channels(1):this_recording_channels(2);
    
    %%% One for indexing SNR of denoised reref data (no NaNs);
    this_recording_rerefs = this_recording_channels(this_recording_channels ~= this_recording.reref_ch);
    
    subject = this_recording.subject{:};
    folder  = this_recording.folder{:};
    bank    = this_recording.bank{:};
    
    subject_ID   = this_recording.subject_ID;
    session_ID   = this_recording.session_ID;
    electrode_ID = this_recording.electrode_ID;
    task_ID      = this_recording.task_ID;
    
    %%% Signal to noise ratios of every channel in the electrode (bank)
    SNR_split = this_recording.SNR_split{:};
    SNR_reref = this_recording.SNR_reref{:};
    SNR_clean = this_recording.SNR_clean{:};
    
    %%% Whether digital or analog headstage was used
    digital_headstage = this_recording.digital_headstage;
    
    %%% Load channels from spike data directory
    spike_directory = fullfile(data_directory, subject, folder, 'spike_data', sprintf('Bank%s', bank));
    
    channels_file = fullfile(spike_directory, 'channels.mat');
    neurons_file  = fullfile(spike_directory, 'neurons.mat');
    classes_file  = fullfile(spike_directory, 'classes.mat');
    
    load(channels_file, 'channels');
    load(neurons_file, 'neurons');
    load(classes_file, 'classes');
    
    %%% Avoid cases where no channels of the electrode yielded any spike
    %%% groups or classes (channels would be an empty array)
    if ~isempty(channels)
        
        %%% Convert channels, neurons, and classes structures to table for concatenation and addition of upper level information
        channels = struct2table(channels);
        channel_numbers = channels.channel_number;
        
        channel_IDs = arrayfun(@(x) num2str(x), channel_numbers, 'UniformOutput', false);
        channel_IDs = strcat(subject, '_', channel_IDs);
        channel_IDs = cell2mat(cellfun(@(x) find(strcmp(unique_channels, x), 1, 'first'), channel_IDs, 'UniformOutput', false));
        
        channels.SNR_split = cell2mat(arrayfun(@(x) SNR_split(this_recording_channels == x), channel_numbers, 'UniformOutput', false));
        channels.SNR_reref = cell2mat(arrayfun(@(x) SNR_reref(this_recording_channels == x), channel_numbers, 'UniformOutput', false));
        channels.SNR_clean = cell2mat(arrayfun(@(x) SNR_clean(this_recording_rerefs == x), channel_numbers, 'UniformOutput', false));
        
        channels.digital_headstage = repmat(digital_headstage, height(channels), 1);
        
        channels.channel_ID = channel_IDs;
       
        neurons = struct2table(neurons);
        
        %%% Dev: Need to fix variable for cases where it might be empty so that it does not become a cell array
        if iscell(neurons.peak_width)
            peak_widths = neurons.peak_width;
            peak_widths = cellfun(@(x) sum(x), peak_widths, 'UniformOutput', false); %%% This takes care of empty cells, converting them to 0.
            neurons.peak_width = [peak_widths{:}]';
        end
        
        classes = struct2table(classes);
        
        %%% Add subject, session, electrode, and channel IDs to tables
        %%% Delete subject, folder, and bank columns
        channels.subject_ID = repmat(subject_ID, height(channels), 1);
        channels.subject = [];
        
        channels.session_ID = repmat(session_ID, height(channels), 1);
        channels.folder = [];
        
        channels.electrode_ID = repmat(electrode_ID, height(channels), 1);
        channels.bank = [];
        
        channels.task_ID = repmat(task_ID, height(channels), 1);
        
        neurons.subject_ID = repmat(subject_ID, height(neurons), 1);
        neurons.subject = [];
        
        neurons.session_ID = repmat(session_ID, height(neurons), 1);
        neurons.folder = [];
        
        neurons.electrode_ID = repmat(electrode_ID, height(neurons), 1);
        neurons.bank = [];
        
        neurons.task_ID = repmat(task_ID, height(neurons), 1);
        
        neurons.channel_ID = arrayfun(@(x) channel_IDs(find(channel_numbers == x, 1, 'first')), neurons.channel_number);     
        
        neurons.n_groups  = NaN(height(neurons), 1);
        neurons.SNR_split = NaN(height(neurons), 1);
        neurons.SNR_reref = NaN(height(neurons), 1);
        neurons.SNR_clean = NaN(height(neurons), 1);
        
        neurons.digital_headstage = repmat(digital_headstage, height(neurons), 1);
        
        for jdx = 1:height(neurons)
        
            channel_ID = neurons(jdx, :).channel_ID;
            sign       = neurons(jdx, :).is_neg;
            option     = neurons(jdx, :).option;
            
            has_channel = channels.channel_ID == channel_ID;
            has_sign    = channels.is_neg == sign;
            has_option  = channels.option == option;
            
            this_channel = channels(has_channel & has_sign & has_option, :);
            
            neurons(jdx, :).n_groups  = this_channel.n_groups;
            neurons(jdx, :).SNR_split = this_channel.SNR_split;
            neurons(jdx, :).SNR_reref = this_channel.SNR_reref;
            neurons(jdx, :).SNR_clean = this_channel.SNR_clean;
            
        end
        
        classes.subject_ID = repmat(subject_ID, height(classes), 1);
        classes.subject = [];
        
        classes.session_ID = repmat(session_ID, height(classes), 1);
        classes.folder = [];
        
        classes.electrode_ID = repmat(electrode_ID, height(classes), 1);
        classes.bank = [];
        
        classes.task_ID = repmat(task_ID, height(classes), 1);
        
        classes.channel_ID = arrayfun(@(x) channel_IDs(find(channel_numbers == x, 1, 'first')), classes.channel_number);
        
        %%% Loop through all classes to add respective channel's n_groups and SNRs
        %%% (based on parent neuron) and also add parent neuron's n_classes, is_SUA,
        %%% is_MUA, and is_noise
        classes.n_groups  = NaN(height(classes), 1);
        classes.SNR_split = NaN(height(classes), 1);
        classes.SNR_reref = NaN(height(classes), 1);
        classes.SNR_clean = NaN(height(classes), 1);
        classes.n_classes = NaN(height(classes), 1);
        
        classes.is_SUA    = false(height(classes), 1);
        classes.is_MUA    = false(height(classes), 1);
        classes.is_noise  = false(height(classes), 1);
        
        classes.digital_headstage = repmat(digital_headstage, height(classes), 1);
        
        for jdx = 1:height(classes)
            
            %%% Get neuron channel ID, session ID, sign and option to find corresponding channel structure
            this_channel_ID = classes(jdx, :).channel_ID;
            this_sign       = classes(jdx, :).is_neg;
            this_option     = classes(jdx, :).option;
            this_group      = classes(jdx, :).group;
            
            has_channel = neurons.channel_ID == this_channel_ID;
            has_sign    = neurons.is_neg == this_sign;
            has_option  = neurons.option == this_option;
            has_group   = neurons.group == this_group;
            
            this_neuron = neurons(has_channel & has_sign & has_option & has_group, :);
            
            classes(jdx, :).n_groups  = this_neuron.n_groups;
            classes(jdx, :).SNR_split = this_neuron.SNR_split;
            classes(jdx, :).SNR_reref = this_neuron.SNR_reref;
            classes(jdx, :).SNR_clean = this_neuron.SNR_clean;
            classes(jdx, :).n_classes = this_neuron.n_classes;
            classes(jdx, :).is_SUA    = this_neuron.is_SUA;
            classes(jdx, :).is_MUA    = this_neuron.is_MUA;
            classes(jdx, :).is_noise  = this_neuron.is_noise;
        
        end
        
        all_channels = [all_channels; channels];
        all_neurons  = [all_neurons; neurons];
        all_classes  = [all_classes; classes];
        
    end
    
end

end


function su_vs_noise_analysis(quality_directory, observation_directory, observation_type, observation_variables, observation_titles, observations)

%%% Figure parameters
figure_width = 1920;
figure_height = 1080; %%%Fit to screen

n_rows = size(observation_variables, 1);
n_columns = size(observation_variables, 2);

[axes_positions, rows, columns] = generate_figure_grid(figure_width, figure_height, n_rows, n_columns);

%%% Colors for single units, multi units, and noise units
su_color = [17, 138, 178] / 255;
mu_color = [252, 158, 79] / 255;
noise_color = [222, 108, 131] / 255;

%%% Different groups' distributions plotted along y-axis. Noise:1, MU:2, SU: 3
y_limits = [0, 4];
y_ticks = 1:3;
y_tick_labels = {'N', 'MU', 'SU'};

%%% Directory and plot name
this_plot_directory = fullfile(quality_directory, observation_directory);
if ~isfolder(this_plot_directory)
    mkdir(this_plot_directory);
end

plot_filename = fullfile(this_plot_directory, strcat(datestr(now, 'yyyy-mm-dd'), '_', observation_type, '_sua-vs-noise'));

%%% LME formula to compare a single feature(var) of single units vs noise groups.
formula = 'var ~ is_SUA + (1|subject_ID) + (1|subject_ID:electrode_ID) + (1|task_ID)';

%%% All analyses plotted in a grid
axes_handles = zeros(n_rows * n_columns, 1);

figure_handle = figure('Units', 'pixels', 'Visible', 'off')
figure_handle.Position(3) = figure_width;
figure_handle.Position(4) = figure_height;

for idx = 1:size(axes_positions, 1)

    this_variable   = observation_variables{rows(idx), columns(idx)};
    this_title = observation_titles{rows(idx), columns(idx)};
    
    %%% Separate table with necessary information for LME analysis and plotting different distributions
    these_observations = observations(:, {'subject_ID', 'electrode_ID', 'task_ID', 'is_SUA', 'is_MUA', 'is_noise'});
    these_observations.var = eval(sprintf('observations.%s', this_variable));
    
    %%% Filter out NaN, inf, or complex numbers resulting from calculations of certain quality metrics. Remove strong outliers.
    these_observations = clean_table(these_observations);
    
    %%% Convert subject_ID and electrode_ID to categorical for LME analysis and exclude is_MUA observations from analysis
    these_observations.subject_ID   = categorical(these_observations.subject_ID);
    these_observations.electrode_ID = categorical(these_observations.electrode_ID);
    these_observations.task_ID      = categorical(these_observations.task_ID);
    
    lme_su_vs_noise = fitlme(these_observations(~these_observations.is_MUA, :), formula);
    
    t_statistic = lme_su_vs_noise.Coefficients.tStat(2);
    p_value     = lme_su_vs_noise.Coefficients.pValue(2);
    
    %%% Gather central tendency measures from each group for box plots
    boxplot_limits = [min(these_observations.var), max(these_observations.var)];
    
    single_units = these_observations(these_observations.is_SUA, :).var;
    sua_stats = get_boxplot_stats(single_units, boxplot_limits);
    
    multi_units = these_observations(these_observations.is_MUA, :).var;
    mua_stats = get_boxplot_stats(multi_units, boxplot_limits);
    
    noise_units = these_observations(these_observations.is_noise, :).var;
    noise_stats = get_boxplot_stats(noise_units, boxplot_limits);
    
    %%% Calculate kernel smoothed probability density functions for violin plots
    [f_sua, xi_sua]     = get_violin(single_units);
    [f_mua, xi_mua]     = get_violin(multi_units);
    [f_noise, xi_noise] = get_violin(noise_units);
    
    %%% Fit x axis to distribution of feature values with some space for labels
    [x_limits, x_ticks, x_tick_labels] = format_axis(these_observations.var);
    
    %%% Axes occupy all available space in corresponding grid cell
    axes_handles(idx) = axes('Parent', figure_handle, 'Units', 'pixels', 'Position', axes_positions(idx, :));
    
    hold on
    
    %%% Title (name of feature) and statistics on top
    text(mean(x_limits), 3.85, this_title, 'FontSize', 16, 'HorizontalAlignment', 'center')
    text(mean(x_limits), 3.7, sprintf('t:%.2f, p:%.3f', t_statistic, p_value), 'FontSize', 14, 'HorizontalAlignment', 'center', 'Color', [0.3, 0.3, 0.3])
    
    %%% Violin and box plots for single units
    if sua_stats.n > 0
        fill([xi_sua, fliplr(xi_sua)], [f_sua + 3, fliplr(3 - f_sua)], su_color, 'FaceAlpha', 0.6, 'EdgeColor', 'none')
        plot([sua_stats.low, sua_stats.high], [3, 3], '-', 'LineWidth', 1.5, 'Color', su_color)
        plot([sua_stats.q1, sua_stats.q3], [3, 3], '-', 'LineWidth', 4.5, 'Color', su_color)
        plot([sua_stats.median, sua_stats.median], [2.9, 3.1], '-', 'LineWidth', 3, 'Color', [0, 1, 0])
    end
    
    %%% Violin and box plots for multi units
    if mua_stats.n > 0
        fill([xi_mua, fliplr(xi_mua)], [f_mua + 2, fliplr(2 - f_mua)], mu_color, 'FaceAlpha', 0.6, 'EdgeColor', 'none');
        plot([mua_stats.low, mua_stats.high], [2, 2], '-', 'LineWidth', 1.5, 'Color', mu_color)
        plot([mua_stats.q1, mua_stats.q3], [2, 2], '-', 'LineWidth', 4.5, 'Color', mu_color)
        plot([mua_stats.median, mua_stats.median], [1.9, 2.1], '-', 'LineWidth', 3, 'Color', [0, 1, 0])
    end
    
    %%% Violin and box plots for noise units
    if noise_stats.n > 0
        fill([xi_noise, fliplr(xi_noise)], [f_noise + 1, fliplr(1 - f_noise)], noise_color, 'FaceAlpha', 0.6, 'EdgeColor', 'none');
        plot([noise_stats.low, noise_stats.high], [1, 1], '-', 'LineWidth', 1.5, 'Color', noise_color)
        plot([noise_stats.q1, noise_stats.q3], [1, 1], '-', 'LineWidth', 4.5, 'Color', noise_color)
        plot([noise_stats.median, noise_stats.median], [0.9, 1.1], '-', 'LineWidth', 3, 'Color', [0 1 0])
    end
    
    %%% X tick labels on the bottom
    text(x_ticks, repmat(0.15, 1, length(x_ticks)), x_tick_labels, 'FontSize', 12, 'HorizontalAlignment', 'center')
    
    %%% Y tick labels to the left
    text(repmat(x_limits(1) + (diff(x_limits) * 0.025), 1, length(y_ticks)), y_ticks, y_tick_labels, 'FontSize', 16)
    
    %%% Group sizes to the right
    group_sizes = arrayfun(@(x) sprintf('N:%d', x), [noise_stats.n, mua_stats.n, sua_stats.n], 'UniformOutput', false);
    text(repmat(x_limits(2) - (diff(x_limits) * 0.1), 1, length(y_ticks)), y_ticks, group_sizes, 'FontSize', 10)
    
    hold off
    
    %%% Empty xticklabels() and yticklabels(). Labels plotted inside plot (vs. default outside) with code above
    xlim(x_limits); xticks(x_ticks); xticklabels([]);
    ylim(y_limits); yticks(y_ticks); yticklabels([]);
    
end

print(plot_filename, '-dpng');
print(plot_filename, '-dsvg');

close all

end


function digital_vs_analog_analysis(quality_directory, observation_directory, observation_type, observation_variables, observation_titles, observations)

%%% Figure parameters
figure_width = 1920;
figure_height = 1080; %%%Fit to screen

n_rows = size(observation_variables, 1);
n_columns = size(observation_variables, 2);

[axes_positions, rows, columns] = generate_figure_grid(figure_width, figure_height, n_rows, n_columns);

%%% Colors for units extracted from recordings of digital vs. analog headstages
digital_color = [125, 206, 19] / 255;
analog_color = [101, 84, 175] / 255;

%%% Distributions of each group plotted along y-axis. A:1, D:2.
y_limits = [0, 3];
y_ticks = 1:2;
y_tick_labels = {'A', 'D'};

%%% Directory and plot name
this_plot_directory = fullfile(quality_directory, observation_directory);
if ~isfolder(this_plot_directory)
    mkdir(this_plot_directory);
end

plot_filename = fullfile(this_plot_directory, strcat(datestr(now, 'yyyy-mm-dd'), '_', observation_type, '_digital-vs-analog'));

%%% LME formula to compare a single feature of digital vs. analog groups
%%% Since subjects usually only had recordings with digital headstage or
%%% analog headstage, only including a random effect for the specific task.
formula = 'var ~ digital_headstage + (1|task_ID)';

%%% All analyses plotted in a grid
axes_handles = zeros(n_rows * n_columns, 1);

figure_handle = figure('Units', 'pixels', 'Visible', 'off')
figure_handle.Position(3) = figure_width;
figure_handle.Position(4) = figure_height;

for idx = 1:size(axes_positions, 1)
    
    this_variable = observation_variables{rows(idx), columns(idx)};
    this_title = observation_titles{rows(idx), columns(idx)};
    
    %%% Separate table with necessary information for LME analysis and plotting different distributions
    these_observations = observations(:, {'task_ID', 'digital_headstage'});
    these_observations.var = eval(sprintf('observations.%s', this_variable));
    
    %%% To filter out NaN, inf, or complex numbers resulting from calculations of certain quality metrics. Remove strong outliers.
    these_observations = clean_table(these_observations);
    
    %%% Task ID converted to categorical for LMEM
    these_observations.subject_ID = categorical(these_observations.task_ID);
    
    lme_digital_vs_analog = fitlme(these_observations, formula);
    
    t_statistic = lme_digital_vs_analog.Coefficients.tStat(2);
    p_value     = lme_digital_vs_analog.Coefficients.pValue(2);
    
    %%% Gather central tendency measures from each group for box plots
    boxplot_limits = [min(these_observations.var), max(these_observations.var)];
    
    digital = these_observations(these_observations.digital_headstage, :).var;
    digital_stats = get_boxplot_stats(digital, boxplot_limits);
    
    analog = these_observations(~these_observations.digital_headstage, :).var;
    analog_stats = get_boxplot_stats(analog, boxplot_limits);
    
    %%% Calculate kernel smoothed probability density functions for violin plots
    [f_digital, xi_digital] = get_violin(digital);
    [f_analog, xi_analog] = get_violin(analog);
    
    %%% Fit x axis to distribution of feature values with some space for labels
    [x_limits, x_ticks, x_tick_labels] = format_axis(these_observations.var);
    
    %%% Axes occupy all available space in corresponding grid cell
    axes_handles(idx) = axes('Parent', figure_handle, 'Units', 'pixels', 'Position', axes_positions(idx, :));
    
    hold on
    
    %%% Title (name of feature) and statistics on top
    text(mean(x_limits), 2.85, this_title, 'FontSize', 16, 'HorizontalAlignment', 'center')
    text(mean(x_limits), 2.7, sprintf('t:%.2f, p:%.3f', t_statistic, p_value), 'FontSize', 14, 'HorizontalAlignment', 'center', 'Color', [0.3, 0.3, 0.3])
    
    %%% Violin and box plots for digital group
    if digital_stats.n > 0
        fill([xi_digital, fliplr(xi_digital)], [f_digital + 2, fliplr(2 - f_digital)], digital_color, 'FaceAlpha', 0.6, 'EdgeColor', 'none');
        plot([digital_stats.low, digital_stats.high], [2, 2], '-', 'LineWidth', 1.5, 'Color', digital_color)
        plot([digital_stats.q1, digital_stats.q3], [2, 2], '-', 'LineWidth', 4.5, 'Color', digital_color);
        plot([digital_stats.median digital_stats.median], [1.9, 2.1], '-', 'LineWidth', 3, 'Color', [1, 0.5, 0])
    end
    
    %%% Violin and box plots for analog group
    if analog_stats.n > 0
        fill([xi_analog, fliplr(xi_analog)], [f_analog + 1, fliplr(1 - f_analog)], analog_color, 'FaceAlpha', 0.6, 'EdgeColor', 'none');
        plot([analog_stats.low, analog_stats.high], [1, 1], '-', 'LineWidth', 1.5, 'Color', analog_color)
        plot([analog_stats.q1, analog_stats.q3], [1, 1], '-', 'LineWidth', 4.5, 'Color', analog_color);
        plot([analog_stats.median, analog_stats.median], [0.9, 1.1], '-', 'LineWidth', 3, 'Color', [1, 0.5, 0])
    end
    
    %%% X tick labels on the bottom
    text(x_ticks, repmat(0.15, 1, length(x_ticks)), x_tick_labels, 'FontSize', 12, 'HorizontalAlignment', 'center')
    
    %%% Y tick labels to the left
    text(repmat(x_limits(1)+(diff(x_limits) * 0.025), 1, length(y_ticks)), y_ticks, y_tick_labels, 'FontSize', 16)
    
    %%% Group sizes to the right
    group_sizes = arrayfun(@(x) sprintf('N:%d', x), [analog_stats.n, digital_stats.n], 'UniformOutput', false);
    text(repmat(x_limits(2) - (diff(x_limits) * 0.1), 1, length(y_ticks)), y_ticks, group_sizes, 'FontSize', 10)
    
    hold off
    
    %%% Empty xticklabels() and yticklabels(). Labels plotted inside rather than outside of plot in code above
    xlim(x_limits); xticks(x_ticks); xticklabels([]);
    ylim(y_limits); yticks(y_ticks); yticklabels([]);
    
end

print(plot_filename, '-dpng');
print(plot_filename, '-dsvg');

close all

end


function count_quality_correlation(quality_directory, observation_directory, observation_type, observation_y, observation_variables, observation_titles, observations, coloring_scheme)

%%% Figure parameters
figure_width = 1920;
figure_height = 1080; %%%Fit to screen

n_rows = size(observation_variables, 1);
n_columns = size(observation_variables, 2);

[axes_positions, rows, columns] = generate_figure_grid(figure_width, figure_height, n_rows, n_columns);

%%% Different coloring schemes. Random colors for each subject, electrode, or task
switch coloring_scheme
    case 'electrodes'
        unique_electrodes = unique(observations.electrode_ID);
        colors = hsv2rgb([rand(length(unique_electrodes), 1), 0.8 + (0.2 * rand(length(unique_electrodes), 1)), 0.7 + (0.3 * rand(length(unique_electrodes), 1))]);
        observations.color = arrayfun(@(x) colors(unique_electrodes == x, :), observations.electrode_ID, 'UniformOutput', false);

    case 'subjects'
        unique_subjects = unique(observations.subject_ID);
        colors = hsv2rgb([rand(length(unique_subjects), 1), 0.8 + (0.2 * rand(length(unique_subjects), 1)), 0.7 + (0.3 * rand(length(unique_subjects), 1))]);
        observations.color = arrayfun(@(x) colors(unique_subjects == x, :), observations.subject_ID, 'UniformOutput', false);

    case 'tasks'
        unique_tasks = unique(observations.task_ID);
        colors = hsv2rgb([rand(length(unique_tasks), 1), 0.8 + (0.2 * rand(length(unique_tasks), 1)), 0.7 + (0.3 * rand(length(unique_tasks), 1))]);
        observations.color = arrayfun(@(x) colors(unique_tasks == x, :), observations.task_ID, 'UniformOutput', false);
end

%%% Directory and plot name
this_plot_directory = fullfile(quality_directory, observation_directory);
if ~isfolder(this_plot_directory)
    mkdir(this_plot_directory);
end

plot_filename = fullfile(this_plot_directory, strcat(datestr(now, 'yyyy-mm-dd'), '_', observation_type, '_count-quality-correlation_', strrep(observation_y, '_', '-'), '_coloring-', coloring_scheme));

%%% LME formula for correlation of y (counts of SUA, MUA, or noise) with
%%% channel or electrode features, considering variability among subjects,
%%% electodes and tasks
if ismember(observation_type, {'AR', 'SR'})
    formula = 'y ~ var + (1|subject_ID) + (1|subject_ID:electrode_ID)';
else
    formula = 'y ~ var + (1|subject_ID) + (1|subject_ID:electrode_ID) + (1|task_ID)';
end

%%% All correlations plotted in a grid 
axes_handles = zeros(n_rows * n_columns, 1);

figure_handle = figure('Units', 'pixels', 'Visible', 'off')
figure_handle.Position(3) = figure_width;
figure_handle.Position(4) = figure_height;

for idx = 1:size(axes_positions, 1)

    this_variable = observation_variables{rows(idx), columns(idx)};
    this_title = observation_titles{rows(idx), columns(idx)};
    
    these_observations = observations(:, {'subject_ID', 'electrode_ID', 'task_ID', 'color'});
    these_observations.y = eval(sprintf('observations.%s', observation_y));
    these_observations.var = eval(sprintf('observations.%s', this_variable));
    
    %%% To filter out NaN, inf, or complex numbers resulting from calculations of certain quality metrics
    these_observations(isnan(these_observations.y), :) = [];
    these_observations = clean_table(these_observations);
    
    %%% Convert subject, electrode, and task IDs to categorical for LME analysis
    these_observations.subject_ID   = categorical(these_observations.subject_ID);
    these_observations.electrode_ID = categorical(these_observations.electrode_ID);
    these_observations.task_ID      = categorical(these_observations.task_ID);
    
    lme_correlation = fitlme(these_observations, formula);
    
    estimate = lme_correlation.Coefficients.Estimate(2);
    p_value  = lme_correlation.Coefficients.pValue(2);
    
    %%% Gather intercept and slope from LMEM for plotting line of best fit
    beta = fixedEffects(lme_correlation);
    x_range = linspace(min(these_observations.var), max(these_observations.var), 100);
    y_predicted = (beta(2) * x_range) + beta(1);
    
    %%% Adjust x and y axes to distributions of var and y
    [x_limits, x_ticks, x_tick_labels] = format_axis(these_observations.var);
    [y_limits, y_ticks, y_tick_labels] = format_axis(these_observations.y);
    
    %%% Axes occupy all available space in corresponding grid cell
    axes_handles(idx) = axes('Parent', figure_handle, 'Units', 'pixels', 'Position', axes_positions(idx, :));
    
    hold on
    
    %%% Plot title and statistics on top
    text(mean(x_limits), y_limits(2) - (diff(y_limits) * 0.025), this_title, 'FontSize', 16, 'HorizontalAlignment', 'center')
    text(mean(x_limits), y_limits(2) - (diff(y_limits) * 0.075), sprintf('e:%.2f, p:%.3f', estimate, p_value), 'FontSize', 14, 'HorizontalAlignment', 'center', 'Color', [0.3, 0.3, 0.3])
    
    %%% All observations with corresponding color and line of best fit
    scatter(these_observations.var, these_observations.y, repmat(8, height(these_observations), 1), cell2mat(these_observations.color));
    plot(x_range, y_predicted, '-', 'LineWidth', 3, 'Color', [0.4, 0.4, 0.4])
    
    %%% X tick labels on the bottom
    text(x_ticks, repmat(y_limits(1) + (diff(y_limits) * 0.025), 1, length(x_ticks)), x_tick_labels, 'FontSize', 12, 'HorizontalAlignment', 'center');
   
    %%% Y tick labels to the left
    text(repmat(x_limits(1) + (diff(x_limits) * 0.025), 1, length(y_ticks)), y_ticks, y_tick_labels, 'FontSize', 12');
    
    hold off
    
    %%% Empty xticklabels and yticklabels. Labels are plotted inside
    %%% rather than outside axes with code above.
    xlim(x_limits); xticks(x_ticks); xticklabels([]);
    ylim(y_limits); yticks(y_ticks); yticklabels([]);
    
end

print(plot_filename, '-dpng');
print(plot_filename, '-dsvg');

close all

end


function digital_vs_analog_counts(quality_directory, observation_directory, observation_type, observation_variables, observation_titles, observations)

%%% Plot parameters. Single smaller plot for this
figure_width = 900;
figure_height = 600;

%%% Colors for channels digital vs. analog headstages
digital_color = [125, 206, 19] / 255;
analog_color = [101, 84, 175] / 255;

%%% Adjust x axis from 0 to max count number
x_limits = [0 , max([observations.n_sua; observations.n_mua;observations.n_noise]) * 1.1];
x_ticks = 0:(diff(x_limits) / 7):x_limits(2);
x_tick_labels = arrayfun(@(x) num2str(round(x)), x_ticks, 'UniformOutput', false);

%%% Distributions of analog and digital groups plotted along y axis.
%%% Counts of noise around y = 1, MUA around 2, SUA around 3
y_limits = [0, length(observation_variables) + 1];
y_ticks = 0.75:0.5:length(observation_variables) + 0.25;
y_tick_labels = repmat({'A', 'D'}, 1, length(observation_variables));
y_position = [{0.75, 1.25}; {1.75, 2.25}; {2.75, 3.25}];

%%% Directory and plot name
this_plot_directory = fullfile(quality_directory, observation_directory);
if ~isfolder(this_plot_directory)
    mkdir(this_plot_directory);
end
plot_filename = fullfile(this_plot_directory, strcat(datestr(now, 'yyyy-mm-dd'), '_', observation_type, '_digital_analog_counts'));

%%% LME formula for contrasting counts between digital and analog groups
if strcmp(observation_directory, 'channels')
    formula = 'var ~ digital_headstage + (1|electrode_ID) + (1|task_ID)';
else
    formula = 'var ~ digital_headstage + (1|task_ID)';
end

figure_handle = figure('Units', 'pixels', 'Visible', 'off')
figure_handle.Position(3) = figure_width;
figure_handle.Position(4) = figure_height;

hold on

%%% Title on top
this_title = 'group type counts';
text(mean(x_limits), y_limits(2) - 0.05 * (diff(y_limits)), this_title, 'FontSize', 16, 'HorizontalAlignment', 'center')

%%% Pre allocate space for t_statistics and p_values of effects to add labels
t_statistic = zeros(height(observations), 1);
p_value     = zeros(height(observations), 1);

for idx = 1:length(observation_variables)
    this_variable = observation_variables{idx};
    
    %%% Separate table to account for possible NaNs in a given measure
    these_observations = observations(:, {'electrode_ID', 'task_ID', 'digital_headstage'});
    these_observations.var = eval(sprintf('observations.%s', this_variable));
    
    %%% Convert task and electrode IDs to categorical for LME analysis
    these_observations(isnan(these_observations.var), :) = [];
    
    these_observations.electrode_ID = categorical(these_observations.electrode_ID);
    these_observations.task_ID      = categorical(these_observations.task_ID);
    
    lme_digital_vs_analog = fitlme(these_observations, formula);
    
    t_statistic(idx) = lme_digital_vs_analog.Coefficients.tStat(2);
    p_value(idx)     = lme_digital_vs_analog.Coefficients.pValue(2);
    
    %%% Gather central tendency measures for digital and analog groups for
    %%% box plots. Not using ksdensity in this case. Plotting observations
    %%% in a y normally distributed around corresponding y_position
    boxplot_limits = [min(these_observations.var), max(these_observations.var)];
    
    digital = these_observations(these_observations.digital_headstage, :);
    digital.y = (randn(height(digital), 1) * 0.06) + y_position{idx, 2};
    digital_stats = get_boxplot_stats(digital.var, boxplot_limits);
    
    analog = these_observations(~these_observations.digital_headstage, :);
    analog.y = (randn(height(analog), 1) * 0.06) + y_position{idx, 1};
    analog_stats = get_boxplot_stats(analog.var, boxplot_limits);
    
    %%% Boxplot for digital
    if digital_stats.n > 0
        scatter(digital.var, digital.y, repmat(8, height(digital), 1), digital_color)
        plot([digital_stats.low, digital_stats.high], [y_position{idx, 2}, y_position{idx, 2}], '-', 'LineWidth', 1.5, 'Color', [0.4, 0.4, 0.4])
        plot([digital_stats.q1, digital_stats.q3], [y_position{idx, 2}, y_position{idx, 2}], '-', 'LineWidth', 4.5, 'Color', [0.4, 0.4, 0.4]);
        plot([digital_stats.median, digital_stats.median], [y_position{idx, 2} - 0.05, y_position{idx, 2} + 0.05], '-w', 'LineWidth', 2)
        plot(digital_stats.median, 2, '.', 'Color', [1, 0, 1])
    end
    text(x_limits(2) - (diff(x_limits) * 0.1), y_position{idx, 2}, sprintf('n = %d', digital_stats.n), 'FontSize', 12);

    %%% Boxplot for analog
    if analog_stats.n > 0
        scatter(analog.var, analog.y, repmat(8, height(analog), 1), analog_color)
        plot([analog_stats.low, analog_stats.high], [y_position{idx, 1}, y_position{idx, 1}], '-', 'LineWidth', 1.5, 'Color', [0.4, 0.4, 0.4])
        plot([analog_stats.q1, analog_stats.q3], [y_position{idx, 1}, y_position{idx, 1}], '-', 'LineWidth', 4.5, 'Color', [0.4, 0.4, 0.4]);
        plot([analog_stats.median, analog_stats.median], [y_position{idx, 1} - 0.05, y_position{idx, 1} + 0.05], '-w', 'LineWidth', 2)
        plot(analog_stats.median, 1, '.', 'Color', [1, 0, 1])
    end
    text(x_limits(2) - (diff(x_limits) * 0.1), y_position{idx, 1}, sprintf('n = %d', analog_stats.n), 'FontSize', 12);
    
end

%%% X tick labels on the bottom
text(x_ticks, repmat(0.15, 1, length(x_ticks)), x_tick_labels, 'FontSize', 12, 'HorizontalAlignment', 'center');

%%% Y tick labels for analog and digital groups
text(repmat(x_limits(1) + (diff(x_limits) * 0.025), 1, length(y_ticks)), y_ticks, y_tick_labels, 'FontSize', 12);

%%% Extra y tick labels for measurement title and statistics placed in between analog and digital groups
text(repmat(x_limits(1) + (diff(x_limits) * 0.05), 1, length(observation_titles)), 1:length(observation_titles), observation_titles, 'FontSize', 14);
stats = arrayfun(@(x) sprintf('t:%.1f, p:%.2f', t_statistic(x), p_value(x)), 1:length(observation_variables), 'UniformOutput', false);
text(repmat(x_limits(2) - (diff(x_limits) * -0.1), 1, length(stats)), 1:length(stats), stats, 'FontSize', 12);

hold off

%%% Empty xticklabels and yticklabels. Plotted inside vs outside with code above
xlim(x_limits); xticks(x_ticks); xticklabels([]);
ylim(x_limits); yticks(x_ticks); yticklabels([]);

print(plot_filename, '-dpng');
print(plot_filename, '-dsvg');

close all

end


function electrode_quality_plots(quality_directory, observation_directory, observation_type, observation_variables, observation_titles, observations, coloring_scheme)

%%% Plot parameters
figure_width = 1920;
figure_height = 1080; %%%Fit to screen

n_rows = length(observation_variables);
subplot_height = floor(figure_height / n_rows);
subplot_width = figure_width;

%%% Coloring schemes. Random colors for each electrode or subjects, or analog vs. digital colors
switch coloring_scheme
    case 'electrodes'
        unique_electrodes = unique(observations.electrode_ID);
        colors = hsv2rgb([rand(length(unique_electrodes), 1), 0.8 + (0.2 * rand(length(unique_electrodes), 1)), 0.7 + (0.3 * rand(length(unique_electrodes), 1))]);
        observations.color = arrayfun(@(x) colors(unique_electrodes == x, :), observations.electrode_ID, 'UniformOutput', false);

    case 'subjects'
        unique_subjects = unique(observations.subject_ID);
        colors = hsv2rgb([rand(length(unique_subjects), 1), 0.8 + (0.2 * rand(length(unique_subjects), 1)), 0.7 + (0.3 * rand(length(unique_subjects), 1))]);
        observations.color = arrayfun(@(x) colors(unique_subjects == x, :), observations.subject_ID, 'UniformOutput', false);

    case 'digital'
        values = [0, 1];
        colors = [101, 84, 175; 125, 206, 19] / 255;
        observations.color = arrayfun(@(x) colors(values == x, :), observations.digital_headstage, 'UniformOutput', false);
end

%%% Directory and plot name
this_plot_directory = fullfile(quality_directory, observation_directory);
if ~isfolder(this_plot_directory)
    mkdir(this_plot_directory);
end

plot_filename = fullfile(this_plot_directory, strcat(datestr(now, 'yyyy-mm-dd'), '_', observation_type, '_electrode-quality_coloring-', coloring_scheme));

%%% All measurements plotted in an n_rows x 1 grid
axes_handles = zeros(n_rows, 1);

figure_handle = figure('Units', 'pixels', 'Visible', 'off')
figure_handle.Position(3) = figure_width;
figure_handle.Position(4) = figure_height;

for idx = 1:n_rows

    this_variable = observation_variables{idx};
    this_title = observation_titles{idx};
    
    these_observations = observations(:, {'color'});
    these_observations.var = eval(sprintf('observations.%s', this_variable));
    
    if ismember(this_variable, {'SNR_split', 'SNR_clean', 'SNR_reref', 'isi_SNR'})
        these_observations.var = cell2mat(cellfun(@(x) mean(x, 'omitnan'), these_observations.var, 'UniformOutput', false));
    end
    
    these_observations = clean_table(these_observations);
        
    x_limits = [-0.1 * height(these_observations), height(these_observations)];
    x_ticks = [1, 50:50:height(these_observations)];
    x_tick_labels = arrayfun(@(x) num2str(x), x_ticks, 'UniformOutput', false);
        
    [y_limits, y_ticks, y_tick_labels] = format_axis(these_observations.var);
    
    axes_handles(idx) = axes('Parent', figure_handle, 'Units', 'pixels', 'Position', [0, figure_handle.Position(4) - (idx * subplot_height), subplot_width, subplot_height]);
    
    hold on
    
    text(mean(x_limits), y_limits(2) - (diff(y_limits)*0.05), this_title, 'FontSize', 16, 'HorizontalAlignment', 'center')
    
    for jdx = 1:height(these_observations)
        bar(jdx, these_observations.var(jdx), 'FaceColor', these_observations.color{jdx})
    end
    
    plot([0, height(these_observations) + 1], [0, 0], '-k')
    
    text(repmat(-0.05 * height(these_observations), 1, length(y_ticks)), y_ticks, y_tick_labels, 'FontSize', 12);
    
    text(x_ticks, repmat(y_limits(1) + (0.05 * diff(y_limits)), 1, length(x_ticks)), x_tick_labels, 'FontSize', 12);
    
    hold off
    
    xlim(y_limits); xticks(y_ticks); xticklabels([]);
    ylim(y_limits); yticks(y_ticks); yticklabels([]);
    
end

print(plot_filename, '-dpng');
print(plot_filename, '-dsvg');

close all

end


function [axes_positions, rows, columns] = generate_figure_grid(figure_width, figure_height, n_rows, n_columns)

axes_positions = zeros(n_rows * n_columns, 4);

subplot_height = floor(figure_height / n_rows);
subplot_width = floor(figure_width / n_columns);

[rows, columns] = ndgrid(1:n_rows, 1:n_columns);

rows = rows(:);
columns = columns(:);

for idx = 1:(n_rows * n_columns)
    axes_positions(idx, :) = [(columns(idx) - 1) * subplot_width, figure_height - (rows(idx) * subplot_height), subplot_width, subplot_height];
end

end


function clean = clean_table(dirty)

dirty(isnan(dirty.var), :) = [];
dirty(isinf(dirty.var), :) = [];

clean = dirty(imag(dirty.var)==0, :);
clean.var = real(clean.var);

zscore = (clean.var - mean(clean.var)) / std(clean.var);

clean(abs(zscore) > 5, :) = [];

end


function boxplot_stats = get_boxplot_stats(observations, boxplot_limits)

q1 = quantile(observations, 0.25);
q3 = quantile(observations, 0.75);

low = q1 - (1.5 * (q3 - q1));
high = q3 + (1.5 * (q3 - q1));

low = max([low, min(boxplot_limits)]);
high = min([high, max(boxplot_limits)]);

boxplot_stats = table;
boxplot_stats.n      = length(observations);
boxplot_stats.median = median(observations);
boxplot_stats.q1     = q1;
boxplot_stats.q3     = q3;
boxplot_stats.low    = low;
boxplot_stats.high   = high;

end


function [f, xi] = get_violin(observations)

f = [];
xi = [];

n = length(observations);

if n > 0

    %%% Rule of thumb calculation of bandwidth
    sigma = std(observations);
    bandwidth = (4 * (sigma^5) / (3 * n) )^(1/5);

    %%% Limit violin to range of observations
    x_range = linspace(min(observations), max(observations), 200);

    %%% Calculate kernel-smoother probability density function
    [f, xi] = ksdensity(observations, x_range, 'Bandwidth', bandwidth);

    %%% Flank probability density functions with zeros and duplicate xi ends so that violin fill is complete
    f = [0, f, 0];
    xi = [xi(1), xi, xi(end)];

    %%% Normalized based on max value so that each violin max width is 0.9
    f = (f / max(f)) * 0.45;
    
end

end


function [limits, ticks, labels] = format_axis(observations)

if diff([min(observations), max(observations)]) == 0

    limits = [min(observations) - 0.5, min(observations) + 0.5];
    ticks = limits(1):0.5:limits(2);
    
else

    limits = [min(observations), max(observations)];
    difference = diff(limits);
    ticks = limits(1):(difference/6):limits(2);
    limits(1) = limits(1) - (difference * 0.15);
    limits(2) = limits(2) + (difference * 0.15);
    
end

labels = arrayfun(@(x) sprintf('%.2f', x), ticks, 'UniformOutput', false);

end