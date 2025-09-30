%%% Base function of Micros Pipeline, Step 3: Filtering Noise
%%% This function filters line noise from rereferenced signal
%%% using 'zapline' by NoiseTools. Calculates signal to noise
%%% ratio after filtering. Makes plots of sample signal and
%%% power spectral density.

function n03_zapline(varargin)
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

%%% List directories
code_directory = fullfile(root_directory, 'micros_pipeline/base_code'); addpath(genpath(code_directory));
check_directory = fullfile(root_directory, 'micros_pipeline/process_files/n03_clean/plots');
data_directory = fullfile(root_directory, 'micros_database');
split_directory = fullfile(data_directory, subject, folder, 'split');
reref_directory = fullfile(data_directory, subject, folder, 'reref', sprintf('Bank%s', bank));
save_directory = strrep(reref_directory, 'reref', 'clean');
if ~exist(save_directory, 'dir')
    mkdir(save_directory)
end
plot_directory = fullfile(check_directory, subject);
if ~exist(plot_directory, 'dir')
    mkdir(plot_directory)
end

tasks_with_stimulation = {'BN', 'FX', 'ARS'};

%%% Get rerefed files and exclude reref_ch
reref_files = dir(fullfile(reref_directory, 'NS6_*.mat'));
reref_files = {reref_files.name};
reref_files = reref_files(~contains(reref_files, {'reref_ch', 'SNR'}));

[file_length, sampling_rate] = get_file_info(split_directory);

%%% Initialize a matrix to hold data
reref_data = zeros(length(reref_files), file_length);

%%% Load files to fill matrix converting data from int16 to double
for idx = 1:length(reref_files)
    reref_file = reref_files{idx};
    load(fullfile(reref_directory, reref_file), 'data', 'sampling_rate'); %%% Loads variables 'data' and 'sampling_rate' (sampling rate)
    reref_data(idx, :) = double(data(1:end));
    clear data  %%% Keep sampling_rate
end

%%% Define filters for zapline method
line_noise_filter = 60 / sampling_rate;

%%% Define filters for stimulation delivered, based on task and subject

if contains(folder, tasks_with_stimulation)
    if contains(folder, 'BN')
        if contains(folder, 'session0')
            stimulation_filter = 100 / sampling_rate;
        elseif contains(folder, 'session1')
            stimulation_filter = 150 / sampling_rate;
        elseif contains(folder, 'session2')
            if strcmp(subject, 'UT289')
                stimulation_filter = 150 / sampling_rate;
            else
                stimulation_filter = 100 / sampling_rate;
            end
        elseif contains(folder, 'session3')
            stimulation_filter = 150 / sampling_rate;
        end
    else
        stimulation_filter = 150 / sampling_rate;
    end
else
    stimulation_filter = [];
end

%%% Other zapline parameters
n_components_removed = 1;            %%% Number of noise components to remove. Default 1. Maybe more would be too much considering one channel removed for referencing.
zapline_parameters.nfft = 65536;     %%% Size of nfft to calculate power. Default 1024. Selected a bit more than 2 seconds (next power of 2 of 60000 samples).
zapline_parameters.nkeep = [];       %%% Number of digital signal components to keep. [] is all. 
zapline_parameters.niterations = 2;  %%% Number of iterations for smoothing filter. Default 1

%%% Clean data using zapline
clean_data = nt_zapline(reref_data', line_noise_filter, n_components_removed, zapline_parameters);
if ~isempty(stimulation_filter)
    clean_data = nt_zapline(clean_data, stimulation_filter, n_components_removed, zapline_parameters);
end
clean_data = clean_data';

%%% Convert clean_data to int16 before saving
clean_data_int = int16(clean_data);

%%% Save clean data to new file
for idx = 1:length(reref_files)  %should only loop through first 8 rows.
    %%% Clean data saved in the same format in /clean directory.
    reref_file = reref_files{idx};
    clean_file = reref_file;
    file_save = fullfile(save_directory, clean_file);
    data = clean_data_int(idx, :);
    save(file_save, 'data', 'sampling_rate')
end

clear clean_data_int

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
SNR_clean = PSD_su - PSD_noise;
save(fullfile(save_directory, 'SNR_clean.mat'), 'SNR_clean');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Plot signals and power spectra

%%% Colors:
reref_color = [55 81 95]/255;   %%% Rerefed channel, unfiltered
clean_color = [23 190 187]/255; %%% Rerefed channel, filtered

%%% Plot parameters
plot_name = fullfile(plot_directory, sprintf('%s_Bank%s', folder, bank));

n_rows = size(reref_data, 1);

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

reref_data_samples = reref_data(:, random_portion);
clean_data_samples = clean_data(:, random_portion);

%%% Butterworth filter. Just to visualize how it would look if filtered by combinato
min_pass = 300; %Hz
max_pass = 3000; %Hz
order = 3;
[butterworth_b, butterworth_a] = butter(order, [min_pass, max_pass]/(sampling_rate / 2), 'bandpass');

%%% Filter plotted samples
%%% Transpose rows to columns, and back to rows
reref_data_samples = filtfilt(butterworth_b, butterworth_a, reref_data_samples')';
clean_data_samples = filtfilt(butterworth_b, butterworth_a, clean_data_samples')';

%%% Calculate appropriate limits to visualize plots. Tick marks spread from 0 by 50 uV (200 because data is not rescaled). 
min_reref = min(reref_data_samples(:));
max_reref = max(reref_data_samples(:));
max_abs = max([abs(max_reref), abs(min_reref)]);
if max_abs < 200
    ylim_reref = 200;
    yticks_reref = [-200, 0, 200];
    yticklabels_reref = {'', '0uV', '50uV'};
else
    ylim_reref = max_abs;
    yticks_reref = sort(unique([(-1 * ylim_reref), (ylim_reref), 0, 0:200:ylim_reref, 0:-200:(-1 * ylim_reref)]));
    yticklabels_reref = repelem({''}, 1, length(yticks_reref));
    yticklabels_reref{yticks_reref == 0} = '0uV';
    yticklabels_reref{yticks_reref == max_abs} = strcat(num2str(round(max_abs / 4)), 'uV'); 
    yticklabels_reref{yticks_reref == 200} = '50uV';
end
ylim_reref = [(-1.1 * ylim_reref), (1.1 * ylim_reref)];

min_clean = min(clean_data_samples(:));
max_clean = max(clean_data_samples(:));
max_abs = max([abs(max_clean), abs(min_clean)]);
if max_abs < 200
    ylim_clean = 200;
    yticks_clean = [-200, 0, 200];
    yticklabels_clean = {'', '0uV', '50uV'};
else
    ylim_clean = max_abs;
    yticks_clean = sort(unique([(-1*ylim_clean), (ylim_clean), 0, 0:200:ylim_clean, 0:-200:(-1*ylim_clean)]));
    yticklabels_clean = repelem({''}, 1, length(yticks_clean));
    yticklabels_clean{yticks_clean == 0} = '0uV';
    yticklabels_clean{yticks_clean == max_abs} = strcat(num2str(round(max_abs/4)), 'uV'); 
    yticklabels_clean{yticks_clean == 200} = '50uV';
end
ylim_clean = [(-1.1*ylim_clean), (1.1*ylim_clean)];

%%% Parameters for generating power spectra for rerefed and clean data
frequency_resolution = 3;                                              %%% Frequency resolution in Hz for power spectrum
sample_size = 2^nextpow2(ceil(sampling_rate/frequency_resolution));    %%% Calculation of number of samples needed to achieve freq_res with fft
frequency_resolution = sampling_rate/sample_size;                      %%% Actual frequency resolution
min_frequency = 3;                                                     %%% Minimum frequency for power spectrum plot
max_frequency = 606;                                                   %%% Maximum frequency for power spectrum plot
frequencies = 0:frequency_resolution:(sampling_rate / 2);              %%% The frequencies corresponding to fft result, 0 to nyquist freq
min_idx = find(frequencies >= min_freq, 1, 'first');                   %%% Indices to use to only store fft results of interest
max_idx = find(frequencies <= max_freq, 1, 'last');
random_portion = random_portion_idx:random_portion_idx + sample_size - 1;

%%% Samples of rerefed and denoised data to generate PSDs
reref_power_samples = reref_data(:, random_portion);
clean_power_samples = clean_data(:, random_portion);

%%% Transpose rows to columns, and back to rows
reref_power_samples = filtfilt(butterworth_b, butterworth_a, reref_power_samples')';
clean_power_samples = filtfilt(butterworth_b, butterworth_a, clean_power_samples')'; 

%%% Convert to decibels
reref_PSD = 10*log10(generate_PSD(reref_power_samples, sampling_rate, min_idx, max_idx));
clean_PSD = 10*log10(generate_PSD(clean_power_samples, sampling_rate, min_idx, max_idx));

%%% Determine appropriate x and y ticks for PSD plots
xlim_PSD = [min_freq, max_freq]; %%% PSDs plotted in decibels (y) and linear freq (x)
xticks_PSD = [min_freq, 60:60:max_freq]; %%% Each frequency tick centered at 60 Hz and harmonics
xticklabels_PSD = repelem({''}, 1, length(xticks_PSD));
xticklabels_PSD{xticks_PSD==60} = '60Hz'; 
min_PSD = min([min(reref_PSD(~isinf(reref_PSD))), min(clean_PSD(~isinf(clean_PSD)))]);
max_PSD = max([max(reref_PSD(~isinf(reref_PSD))), max(clean_PSD(~isinf(clean_PSD)))]);
if max_PSD < 10 || isempty(max_PSD)
    max_PSD = 10;
end
if min_PSD > 0 || isempty(min_PSD)
    min_PSD = 0;
end  

ylim_PSD = [min_PSD, max_PSD];
yticks_PSD = sort(unique([min_PSD, max_PSD, 0, 0:10:max_PSD, 0:-10:min_PSD]));
ylim_PSD = [min(ylim_PSD)-(diff(ylim_PSD)*0.05), max(ylim_PSD)+(diff(ylim_PSD)*0.05)];
yticklabels_PSD = repelem({''}, 1, length(yticks_PSD));
yticklabels_PSD{yticks_PSD == 10} = '10dB';
yticklabels_PSD{yticks_PSD == 0} = '0dB';
if round(max_PSD) > 10
    yticklabels_PSD{yticks_PSD == max_PSD} = strcat(num2str(max_PSD), 'dB');
end
if round(min_PSD) < 10 && round(min_PSD) ~= 0
    yticklabels_PSD{yticks_PSD == min_PSD} = strcat(num2str(min_PSD), 'dB');
end

%%% Initialize figure to be the size of the screen and initialize matrix to hold axes for 3 subplots per channel (n_rows), the rerefed signal, the denoised signal, and PSDs for both.
figure('Units', 'pixels', 'Visible', 'off')
figure_handle = gcf;
figure_handle.Position(3) = figure_width;
figure_handle.Position(4) = figure_height;
axes_handles = zeros(n_rows, 3);

for idx = 1:n_rows
    %%% One color for rerefed signal and one color for denoised signal
    
    %%% Plot of rerefed signal, 250ms sample
    axes_handles(idx, 1) = axes('Parent', figure_handle, 'Units', 'pixels', 'Position', [0, figure_handle.Position(4) - (idx * subplot_height), subplot_width1, subplot_height]);
    
    plot(reref_data_samples(idx, :), 'Color', reref_color, 'Line_width', 1.25)
    
    xlim(xlim1); xticks(xticks1); xticklabels([]);
    ylim(ylim_reref); yticks(yticks_reref); yticklabels([]);
    
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
    
    %%% Plot of denoised signal, same 250ms sample
    axes_handles(idx, 2) = axes('Parent', figure_handle, 'Units', 'pixels', 'Position', [subplot_width1, figure_handle.Position(4) - (idx * subplot_height), subplot_width1, subplot_height]);
    
    plot(clean_data_samples(idx, :), 'Color', clean_color, 'Line_width', 1.25)
    
    xlim(xlim1); xticks(xticks1); xticklabels([]);
    ylim(ylim_clean); yticks(yticks_clean); yticklabels([])
    
    if idx == n_rows
        for jdx = 1:length(xticks1)
            if ~isempty(xticklabels1{jdx})
                text(xticks1(jdx), min(ylim_clean)+(0.05*diff(ylim_clean)), xticklabels1{jdx}, 'Color', [0.7 0 0.7], 'FontSize', 16, 'HorizontalAlignment', 'center');
            end
        end
        for jdx = 1:length(yticks_clean)
            if ~isempty(yticklabels_clean{jdx})
                text(sampling_rate * 0.125, yticks_clean(jdx), yticklabels_clean{jdx}, 'Color', [0.7 0 0.7], 'FontSize', 16, 'HorizontalAlignment', 'center');
            end
        end
    end 
    
    %%% Plot of PSDs for rerefed and denoised signal
    axes_handles(idx, 3) = axes('Parent', figure_handle, 'Units', 'pixels', 'Position', [subplot_width1 * 2, figure_handle.Position(4) - (idx * subplot_height), subplot_width2, subplot_height]);
    
    plot(freqs(min_idx:max_idx), reref_PSD(idx, :), 'Color', reref_color)
    hold on
    plot(freqs(min_idx:max_idx), clean_PSD(idx, :), 'Color', clean_color)
    hold off
    
    xlim(xlim_PSD); xticks(xticks_PSD); xticklabels([]); 
    ylim(ylim_PSD); yticks(yticks_PSD); yticklabels([])
    
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
PSD = zeros(size(data_samples, 1), max_idx-min_idx+1);
N = size(data_samples, 2);
for idx = 1:size(data_samples, 1)
    sample = data_samples(idx, :);
    X = fft(sample); %%%Calculation of fourier coefficients
    Pxx = (abs(X(1:N/2+1)).^2) * 2 / (sampling_rate * N); %%%PSD calculation of real signal
    Pxx(1) = Pxx(1) / 2; %%Correction of 0Hz DC value
    PSD(idx, :) = Pxx(min_idx:max_idx); %%%Only keeping the freq range to be plotted
end
end


function [file_length, sampling_rate] = get_file_info(split_directory)
measuresFile = fullfile(split_directory, 'file_length.mat');
measures = load(measuresFile);
file_length = measures.file_length(1);
sampling_rate = measures.sampling_rate(1);
end