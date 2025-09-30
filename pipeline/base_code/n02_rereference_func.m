%%% Base function of Micros Pipeline, Step 2: Rereferencing Recording
%%% Based on selected mode of referencing, this code will rereference
%%% all channel data of a single microelectrode and save copies.
%%% Additionally calculates signal to noise ratios before and after
%%% referencing. Makes plots of sample signal and power spectral density.

function n02_rereference_func(varargin)
if isempty(varargin)                                               %%% To run manually, edit values below
    root_directory = '/path/to/micros_pipeline/parent_directory';  %%% Root directory with pipeline and database folders
    subject = 'SC000';                                             %%% Subject code in SC000 format
    folder = 'yyyy-mm-dd_task-code_part1';                         %%% Folder in yyyy-mm-dd_task-code format, or yyyy-mm-dd_task-code_part1 format if more than 1 part.
    bank = 'A'                                                     %%% Letter of microelectrode bank ('A', 'B', 'C', 'D')
    channels = {[1, 8]};                                           %%% [min channel number of bank, max channel number of bank]
    has_sync_pulses = true;                                        %%% Whether there was sync channel after splitting recording (true or false).
    sync_channel = 129;                                            %%% Channel number containing sync pulse data.
    referencing_mode = 'signal_power';                             %%% Mode of referencing: 
                                                                   %%% 'common_average' (common average referencing), 
                                                                   %%% 'signal_power' (uses channel with lowest signal power as reference),
                                                                   %%% 'none' (only highpasses data over 0.5Hz and demeans it)
else                                                               %%% Otherwise, this is the order they should be entered into function, following above format
    root_directory = varargin{1};
    subject = varargin{2};
    folder = varargin{3};
    bank = varargin{4};
    channels = varargin{5};
    has_sync_pulses = varargin{6};
    sync_channel = varargin{7};
    referencing_mode = varargin{8};
end

%%% List directories and add code directory to path
code_directory = fullfile(root_directory, 'micros_pipeline/base_code'); addpath(genpath(code_directory));
check_directory = fullfile(root_directory, 'micros_pipeline/process_files/n02_reref/plots');
data_directory = fullfile(root_directory, 'micros_database');
split_directory = fullfile(data_directory, subject, folder, 'split');
save_directory = strrep(split_directory, 'split', 'reref');
save_directory = fullfile(save_directory, sprintf('Bank%s', bank));
if ~exist(save_directory, 'dir')
    mkdir(save_directory)
end
plot_directory = fullfile(check_directory, subject);
if ~exist(plot_directory, 'dir')
    mkdir(plot_directory)
end

%%% Get split files and exclude sync channel, jacksheet, params and time stamps file
split_files = dir(fullfile(split_directory, 'NS6_*.mat'));
split_files = {split_files.name};
if has_sync_pulses
    split_files = split_files(~contains(split_files, {sprintf('_%03d', sync_channel), 'time_stamps', 'jacksheet', 'params', 'length'}));
else
    split_files = split_files(~contains(split_files, {'time_stamps', 'jacksheet', 'params', 'length'}));
end

%%% Only include channels belonging to the microelectrode bank    
other_microelectrode_bank = false(length(split_files), 1);
for idx = 1:length(split_files)
    split_file = split_files{idx};
    channel_number = str2double(split_file(end - 6:end - 4));
    if channel_number < min(channels{:}) || channel_number > max(channels{:})
        other_microelectrode_bank(idx) = true;
    end
end
split_files = split_files(~other_microelectrode_bank);

%%% Get file length and sampling rate
[file_length, sampling_rate] = get_file_info(split_directory);

%%% Initialize matrix to hold data of every channel in microelectrode
split_data = zeros(length(split_files), file_length);

%%% Load files to fill matrix converting data from int16 to double
for idx = 1:length(split_files)
    split_file = split_files{idx};
    load(fullfile(split_directory, split_file), 'data', 'sampling_rate');
    split_data(idx, :) = double(data(1:end));
    clear data  %%% Keep sampling rate
end

%%% For each raw file, will apply highpass filter over 0.5 Hz to minimize
%%% drift and will remove DC offset by demeaning the signal, to achieve
%%% better referencing.

%%% Butterworth filter to highpass data
highpass_frequency = 0.5; %Hz
order = 2;
[butterworth_b, butterworth_a] = butter(order, highpass_freq/(sampling_rate/2), 'high');
split_data = filtfilt(butterworth_b, butterworth_a, split_data')'; %%% Transpose rows to columns, and back to rows

%%% Demean the data
split_data = split_data - mean(split_data, 2);

%%% List of parameters to measure PSD within frequency bands used for single unit detection (signal power)
frequency_resolution_su = 50;                                               %%% Frequency resolution
min_frequency_su = 350;                                                     %%% Minimum frequency for spike detection in combinato
max_frequency_su = 2950;                                                    %%% Maximum frequency for spike detection in combinato
sample_size_su = 2^nextpow2(ceil(sampling_rate/frequency_resolution_su));   %%% Calculation of number of samples needed to achieve freq_res with fft in powers of 2.
frequency_resolution_su = sampling_rate/sample_size_su;                     %%% Actual frequency of resolution. Should be equal or lower (better) than set.
frequencies_su = 0:frequency_resolution_su:(sampling_rate/2);               %%% The frequencies corresponding to fft result, 0 to nyquist freq
min_idx_su = find(frequencies_su >= min_frequency_su, 1, 'first');          %%% Indices to use to only store fft results of interest
max_idx_su = find(frequencies_su <= max_frequency_su, 1, 'last');
    
%%% Calculate power spectral densities (PSDs) for every channel
start_idx = floor(size(split_data, 2) * 0.5);                               %%% Getting power spectrum of halfway through file
samples_su = double(split_data(:, start_idx:(start_idx+sample_size_su-1)));
PSD_su = generate_PSD(samples_su, sampling_rate, min_idx_su, max_idx_su);
PSD_su = mean(PSD_su, 2);
PSD_su = 10*log10(PSD_su);

%%% List of parameters to measure PSD in frequency band of line noise
frequency_resolution_noise = 1;                                                    %%% Frequency resolution
min_frequency_noise = 55;                                                          %%% Minimum frequency for broadband linenoise
max_frequency_noise = 65;                                                          %%% Maximum frequency for broadband linenoise
sample_size_noise = 2^nextpow2(ceil(sampling_rate/frequency_resolution_noise));    %%% Calculation of number of samples needed to achieve freq_res with fft in powers of 2.
frequency_resolution_noise = sampling_rate/sample_size_noise;                      %%% Actual frequency of resolution. Should be equal or lower (better) than set.
frequencies_noise = 0:frequency_resolution_noise:(sampling_rate/2);                %%% The frequencies corresponding to fft renoiselt, 0 to nyquist freq
min_idx_noise = find(frequencies_noise >= min_frequency_noise, 1, 'first');        %%% Indices to use to only store fft renoiselts of interest
max_idx_noise = find(frequencies_noise <= max_frequency_noise, 1, 'last');
    
%%% Calculate PSDs for every channel
samples_noise = double(split_data(:, start_idx:(start_idx+sample_size_noise-1)));
PSD_noise = generate_PSD(samples_noise, sampling_rate, min_idx_noise, max_idx_noise);
PSD_noise = mean(PSD_noise, 2);
PSD_noise = 10*log10(PSD_noise);

%%% Signal-to-noise ratio (SNR) is difference of decibels between signal PSD and noise PSD
SNR_split = PSD_su - PSD_noise;
save(fullfile(save_directory, 'SNR_split.mat'), 'SNR_split');

switch referencing_mode
    case 'signal_power'
        %%% Identify channel with lowest SU power (that seemed to work better than reref to channel with greatest amount of line noise, or lower SNR, or CAR
        %%% because line noise dominates the raw signal and one channel could have significantly more than others)
        %%% The result is more channels with lower magnitude line noise
        
        [~, reref_idx] = min(PSD_su);
        
        %%% Rereference the data based on reref_idx
        reference = split_data(reref_idx, :);
        
        %%% Subtract reref_ch from all data. The row for reref_ch should now be
        %%% zeros.
        rereferenced_data = split_data - reference;
        
    case 'common_average'
        reref_idx = 9;                  %%% Store in 9th row for visualization
        
        CAR = mean(split_data, 1);      %%% Get common average reference
        
        split_data(reref_idx, :) = CAR;
        
        %%% Subtract CAR from all data
        rereferenced_data = split_data - CAR;
        
        %%% Create text file to specify it was CAR and code channel number to 000
        reref_file = fullfile(save_directory, 'reref_ch_000_common_average');    
        file_id = fopen(reref_file, 'w');fclose(file_id);
 
    case 'none'
        rereferenced_data = split_data;
        reref_idx = 0;
        
        %%% Create text file to specify it was none and code channel number to 000
        reref_file = fullfile(save_directory, 'reref_ch_000_none');              
        file_id = fopen(reref_file, 'w'); fclose(file_id);

end

%%% Convert rereferenced_data to int16 before saving
rereferenced_data_int = int16(rereferenced_data);

%%% Save rereferenced data to new file
for idx = 1:length(split_files)
    
    split_file = split_files{idx};
    
    if idx == reref_idx
    
        %%% Create text file to identify channel number used as reference
        reref_channel_number = regexp(split_file, '[0-9][0-9][0-9]', 'match'); %%%Get 3 digits of channel number
        reref_channel_number = reref_channel_number{1};
        reref_file = fullfile(save_directory, sprintf('reref_ch_%s_%s', reref_channel_number, referencing_mode));           
        file_id = fopen(reref_file, 'w'); fclose(file_id);
        
    else
    
        %%% Rerefed data saved in the same format in /reref directory.
        reref_file = split_file;
        file_save = fullfile(save_directory, reref_file); 
        data = rereferenced_data_int(idx, :);
        save(file_save, 'data', 'sampling_rate')
        
    end
    
end

clear rereferenced_data_int

%%% Measure PSDs and SNR for rereferenced_data
samples_su = double(rereferenced_data(:, start_idx:(start_idx+sample_size_su-1)));
PSD_su = generate_PSD(samples_su, sampling_rate, min_idx_su, max_idx_su);
PSD_su = mean(PSD_su, 2);
PSD_su = 10*log10(PSD_su);
    
samples_noise = double(rereferenced_data(:, start_idx:(start_idx+sample_size_noise-1)));
PSD_noise = generate_PSD(samples_noise, sampling_rate, min_idx_noise, max_idx_noise);
PSD_noise = mean(PSD_noise, 2);
PSD_noise = 10*log10(PSD_noise);

%%% SNR is difference of decibels between signal (su) PSD and noise PSD
SNR_reref = PSD_su - PSD_noise;
save(fullfile(save_directory, 'SNR_reref.mat'), 'SNR_reref');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Plot signals and power spectra

%%% Colors:
ref_raw_color = [242 76 0]/255;      %%% reference, raw
ref_cln_color = [252 122 30]/255;    %%% reference, rereferenced
raw_color = [55 81 95]/255;          %%% split channel, raw
cln_color = [23 190 187]/255;        %%% split channel, rereferenced

%%% Plot parameters
plot_name = fullfile(plot_directory, sprintf('%s_Bank%s_%s', folder, bank, mode));

n_rows = size(split_data, 1);

figure_width = 4480;
figure_height = 2520;                                   %%% Fit to screen
subplot_height = round(figure_height / n_rows);
subplot_width1 = (figure_width - subplot_height) / 2;   %%% For sample signal
subplot_width2 = subplot_height;                        %%% For PSD

%%% Plotting quarter second of sample signal for each channel
xlim1 = [0, (sampling_rate / 4)];
xticks1 = 0:(sampling_rate / 20):(sampling_rate / 4);
xticklabels1 = {'0ms', '50ms', '100ms', '150ms', '200ms', '250ms'};

%%% For plotting EEG samples
random_portion_idx = round(size(split_data, 2) * (0.3 + (0.4 * rand)));            %%% Signal randomly sampled somewhere between 30-70% of the recording.
random_portion = random_portion_idx:random_portion_idx + (sampling_rate / 4) - 1;

split_data_samples = split_data(:, random_portion);
rereferenced_data_samples = rereferenced_data(:, random_portion);

%%% Butterworth filter to visualize how signal would look if filtered by combinato
min_pass = 300; %Hz
max_pass = 3000; %Hz
order = 3;
[butterworth_b, butterworth_a] = butter(order, [min_pass, max_pass]/(sampling_rate / 2), 'bandpass');

%%% Filter plotted samples
split_data_samples = filtfilt(butterworth_b, butterworth_a, split_data_samples')'; %%% Transpose rows to columns, and back to rows
rereferenced_data_samples = filtfilt(butterworth_b, butterworth_a, rereferenced_data_samples')';

%%% Calculate appropriate limits to visualize plots. Tick marks spread from 0
%%% by 50 uV (200 because data is not rescaled). For raw data tick marks are
%%% spread from 0 by 500 uV (2000) because data can be that large. However if
%%% the raw data file is flatter, yticks spread by 50.

min_split = min(split_data_samples(:));
max_split = max(split_data_samples(:));
max_abs = max([abs(max_split), abs(min_split)]);

if max_abs < 200
    ylim_split = 200;
    yticks_split = [-200, 0, 200];
    yticklabels_split = {'', '0uV', '50uV'};
else
    ylim_split = max_abs;
    yticks_split = sort(unique([(-1 * ylim_split), (ylim_split), 0, 0:2000:ylim_split, 0:-2000:(-1 * ylim_split)]));
    yticklabels_split = repelem({''}, 1, length(yticks_split));
    yticklabels_split{yticks_split == 0} = '0uV';
    yticklabels_split{yticks_split == max_abs} = strcat(num2str(round(max_abs/4)), 'uV');
    if any(yticks_split == 2000)
        yticklabels_split{yticks_split == 2000} = '500uV';
    end
end
ylim_split = [(-1.1 * ylim_split), (1.1 * ylim_split)];

min_reref = min(rereferenced_data_samples(:));
max_reref = max(rereferenced_data_samples(:));
max_abs = max([abs(max_reref), abs(min_reref)]);

if max_abs < 200
    ylim_reref = 200;
    yticks_reref = [-200, 0, 200];
    yticklabels_reref = {'', '0uV', '50uV'};
else
    ylim_reref = max_abs;
    yticks_reref = sort(unique([(-1*ylim_reref), (ylim_reref), 0, 0:200:ylim_reref, 0:-200:(-1*ylim_reref)]));
    yticklabels_reref = repelem({''}, 1, length(yticks_reref));
    yticklabels_reref{yticks_reref == 0} = '0uV';
    yticklabels_reref{yticks_reref == max_abs} = strcat(num2str(round(max_abs/4)), 'uV'); 
    yticklabels_reref{yticks_reref == 200} = '50uV';
end
ylim_reref = [(-1.1 * ylim_reref), (1.1 * ylim_reref)];

%%% Parameters for generating power spectra for split and reref data
frequency_resolution = 3;                                              %%% Frequency resolution in Hz for power spectrum
sample_size = 2^nextpow2(ceil(sampling_rate/frequency_resolution));    %%% Calculation of number of samples needed to achieve freq_res with fft
frequency_resolution = sampling_rate/sample_size;                      %%% Actual frequency resolution
min_frequency = 3;                                                     %%% Minimum frequency for power spectrum plot
max_frequency = 606;                                                   %%% Maximum frequency for power spectrum plot
frequencies = 0:frequency_resolution:(sampling_rate / 2);              %%% The frequencies corresponding to fft result, 0 to nyquist freq
min_idx = find(frequencies >= min_freq, 1, 'first');                   %%% Indices to use to only store fft results of interest
max_idx = find(frequencies <= max_freq, 1, 'last');
random_portion = random_portion_idx:random_portion_idx + sample_size - 1;

%%% Samples of data to generate PSDs
split_power_samples = split_data(:, random_portion);
reref_power_samples = rereferenced_data(:, random_portion);

%%% Transpose rows to columns, and back to rows
split_power_samples = filtfilt(butterworth_b, butterworth_a, split_power_samples')';
reref_power_samples = filtfilt(butterworth_b, butterworth_a, reref_power_samples')';

%%% Convert to decibels
split_PSD = 10*log10(generate_PSD(split_power_samples, sampling_rate, min_idx, max_idx));
reref_PSD = 10*log10(generate_PSD(reref_power_samples, sampling_rate, min_idx, max_idx));

%%% Determine appropriate x and y ticks for PSD plots
%%% PSDs plotted in decibels (y) and linear freq (x)
%%% Each frequency tick centered at 60 Hz and harmonics

xlim_PSD = [min_freq, max_freq];
xticks_PSD = [min_freq, 60:60:max_freq];
xticklabels_PSD = repelem({''}, 1, length(xticks_PSD));
xticklabels_PSD{xticks_PSD==60} = '60Hz';

min_PSD = min([min(split_PSD(~isinf(split_PSD))), min(reref_PSD(~isinf(reref_PSD)))]);
max_PSD = max([max(split_PSD(~isinf(split_PSD))), max(reref_PSD(~isinf(reref_PSD)))]);
if max_PSD < 10 || isempty(max_PSD)
    max_PSD = 10;
end
if min_PSD > 0 || isempty(min_PSD)
    min_PSD = 0;
end  

ylim_PSD = [min_PSD, max_PSD];
yticks_PSD = sort(unique([min_PSD, max_PSD, 0, 0:10:max_PSD, 0:-10:min_PSD]));
ylim_PSD = [min(ylim_PSD) - (diff(ylim_PSD) * 0.05), max(ylim_PSD) + (diff(ylim_PSD) * 0.05)];
yticklabels_PSD = repelem({''}, 1, length(yticks_PSD));
yticklabels_PSD{yticks_PSD == 10} = '10dB';
yticklabels_PSD{yticks_PSD == 0} = '0dB';

if round(max_PSD) > 10
    yticklabels_PSD{yticks_PSD == max_PSD} = strcat(num2str(round(max_PSD)), 'dB');
end
if round(min_PSD) < 10 && round(min_PSD) ~= 0
    yticklabels_PSD{yticks_PSD == min_PSD} = strcat(num2str(round(min_PSD)), 'dB');
end

%%% Initialize figure to be the size of the screen and initialize matrix to
%%% hold axes for 3 subplots per channel(n_rows), the unrerefed(raw) signal, the rerefed signal, and PSDs for both.

figure('Units', 'pixels', 'Visible', 'off')
figure_handle = gcf;
figure_handle.Position(3) = figure_width;
figure_handle.Position(4) = figure_height;

axes_handles = zeros(n_rows, 3);

for idx = 1:n_rows
    
    %%% Different colors for reref and split channel data, before and after reref
    if idx == reref_idx
        color_raw = ref_raw_color;
        color_clean = ref_cln_color;
    else
        color_raw = raw_color;
        color_clean = cln_color;
    end
    
    %%% Plot of unrerefed (raw) signal, 250ms sample
    axes_handles(idx, 1) = axes('Parent', figure_handle, 'Units', 'pixels', 'Position', [0, figure_handle.Position(4) - (idx * subplot_height), subplot_width1, subplot_height]);
    
    plot(split_data_samples(idx, :), 'Color', color_raw, 'LineWidth', 1.25)
    
    xlim(xlim1)
    ylim(ylim_split)
    xticks(xticks1)
    yticks(yticks_split)
    xticklabels([])
    yticklabels([])
    
    if idx == n_rows
        for jdx = 1:length(xticks1)
            if ~isempty(xticklabels1{jdx})
                text(xticks1(jdx), min(ylim_split)+(0.05*diff(ylim_split)), xticklabels1{jdx}, 'Color', [0.7 0 0.7], 'FontSize', 16, 'HorizontalAlignment', 'center');
            end
        end
        for jdx = 1:length(yticks_split)
            if ~isempty(yticklabels_split{jdx})
                text(sampling_rate*0.125, yticks_split(jdx), yticklabels_split{jdx}, 'Color', [0.7 0 0.7], 'FontSize', 16, 'HorizontalAlignment', 'center');
            end
        end
    end
    
    %%% Plot of rerefed (clean) signal, same 250ms sample
    axes_handles(idx, 2) = axes('Parent', figure_handle, 'Units', 'pixels', 'Position', [subplot_width1, figure_handle.Position(4) - (idx * subplot_height), subplot_width1, subplot_height]);
    
    plot(rereferenced_data_samples(idx, :), 'Color', color_clean, 'LineWidth', 1.25)
    
    xlim(xlim1)
    ylim(ylim_reref)
    xticks(xticks1)
    yticks(yticks_reref)
    xticklabels([])
    yticklabels([])
    
    if idx == n_rows
        for jdx = 1:length(xticks1)
            if ~isempty(xticklabels1{jdx})
                text(xticks1(jdx), min(ylim_reref)+(0.05*diff(ylim_reref)), xticklabels1{jdx}, 'Color', [0.7 0 0.7], 'FontSize', 16, 'HorizontalAlignment', 'center');
            end
        end
        for jdx = 1:length(yticks_reref)
            if ~isempty(yticklabels_reref{jdx})
                text(sampling_rate*0.125, yticks_reref(jdx), yticklabels_reref{jdx}, 'Color', [0.7 0 0.7], 'FontSize', 16, 'HorizontalAlignment', 'center');
            end
        end
    end
    
    %%% Plot of PSDs of rerefed (clean) and unrerefed (raw) signal
    axes_handles(idx, 3) = axes('Parent', figure_handle, 'Units', 'pixels', 'Position', [subplot_width1 * 2, figure_handle.Position(4) - (idx * subplot_height), subplot_width2, subplot_height]);
    
    plot(freqs(min_idx:max_idx), split_PSD(idx, :), 'Color', color_raw)
    hold on
    plot(freqs(min_idx:max_idx), reref_PSD(idx, :), 'Color', color_clean)
    hold off
    
    xlim(xlim_PSD)
    ylim(ylim_PSD)
    xticks(xticks_PSD)
    yticks(yticks_PSD)
    xticklabels([])
    yticklabels([])
    
    if idx == n_rows
        for jdx = 1:length(xticks_PSD)
            if ~isempty(xticklabels_PSD{jdx})
                text(xticks_PSD(jdx), min(ylim_PSD)+(0.05*diff(ylim_PSD)), xticklabels_PSD{jdx}, 'Color', [0.7 0 0.7], 'FontSize', 16, 'HorizontalAlignment', 'center');
            end
        end
        for jdx = 1:length(yticks_PSD)
            if ~isempty(yticklabels_PSD{jdx})
                text(30, yticks_PSD(jdx), yticklabels_PSD{jdx}, 'Color', [0.7 0 0.7], 'FontSize', 16, 'HorizontalAlignment', 'center');
            end
        end
    end
    
end

print(plot_name, '-dpng');
print(plot_name, '-dsvg');
close all

end


function PSD = generate_PSD(data_samples, sampling_rate, min_idx, max_idx)
PSD = zeros(size(data_samples, 1), max_idx - min_idx + 1);
N = size(data_samples, 2);
for idx = 1:size(data_samples, 1)
    sample = data_samples(idx, :);
    X = fft(sample); %%% Calculation of fourier coefficients
    Pxx = (abs(X(1:N/2+1)).^2) * 2 / (sampling_rate * N); %%% PSD calculation of real signal
    Pxx(1) = Pxx(1) / 2; %%% Correction of 0Hz DC value
    PSD(idx, :) = Pxx(min_idx:max_idx); %%% Only keeping the freq range to be plotted
end
end


function [file_length, sampling_rate] = get_file_info(split_directory)
measures_file = fullfile(split_directory, 'file_length.mat');
measures = load(measures_file);
file_length = measures.file_length(1);
sampling_rate = measures.sampling_rate(1);
end