%%% Function to translate .h5 results from combinato clustering.
%%% Uses this information to make a structure for the session containing
%%% quality metrics at the bank, channel, group, and class levels. This
%%% information could be used to automatically classify groups,
%%% refine groups, and refine classes.

function n08_spike_data_func(varargin)
if isempty(varargin)                                               %%% To run manually edit values below
    root_directory = '/path/to/micros_pipeline/parent_directory';  %%% Root directory with pipeline and database folders
    subject = 'SC000';                                             %%% Subject code in SC000 format
    folder = 'yyyy-mm-dd_task-code_part1';                         %%% Folder in yyyy-mm-dd_task-code format, or yyyy-mm-dd_task-code_part1 format if more than 1 part. 
    bank = 'A';                                                    %%% Recording hardware bank character ('A', 'B', 'C', 'D')
    clustering_mode = 'standard';                                  %%% Clustering mode used:
                                                                   %%% 'standard' only info of default (options00)
                                                                   %%% 'grid-search' default + options01-16
else                                                               %%% Otherwise this is the order they should be entered into function, following above format
    root_directory = varargin{1};
    subject = varargin{2};
    folder = varargin{3};
    bank = varargin{4};
    clustering_mode = varargin{5};
end

warning off

%%% List directories
data_directory = fullfile(root_directory, 'micros_database');
split_directory = fullfile(data_directory, subject, folder, 'split');
rescale_directory = fullfile(data_directory, subject, folder, 'rescaled', sprintf('Bank%s', bank));
combinato_directory = fullfile(data_directory, subject, folder, 'combinato_files', sprintf('Bank%s', bank));
save_directory = strrep(combinato_directory, 'combinato_files', 'spike_data');
if ~exist(save_directory, 'dir')
    mkdir(save_directory);
end

%%% Get file info (file_length, sampling_rate, noise_stds) from split_directory and rescale_directory 
file_info = load(fullfile(split_directory, 'file_length.mat'));
sampling_rate = file_info.sampling_rate;
file_length = file_info.file_length;
noise_info = load(fullfile(rescale_directory, 'noise_info.mat'), 'noise_info');
noise_standard_deviations = noise_info.noise_standard_deviation;

%%% Loop through each channel/option/sign combination and convert .h5 file with results of sorting to
%%% matlab variables from which to extract spike data
%%% The results for each channel are stored in two .h5 files:
%%% 
%%% -data_(channel).h5:
%%%   -Single file with negative and positive spikes
%%%   -3 datasets grouped by /neg or /pos depending on polarity of spikes:
%%%
%%%     -/spikes              :       60 x n_spikes   : (single) 60 sample waveforms (2ms)
%%%     -/times               : n_spikes x 1          : (double) the times of n detected spikes
%%%     -/artifacts           : n_spikes x 1          : (int8) logical values pointing to artifactual spikes.
%%% 
%%% -sort_cat.h5:
%%%   -Separate files for clusterings of each channel and negative and positive spikes
%%%   -Does not match size of above array, some neurons skipped. 10 datasets:
%%%   
%%%     -types                :        2 x n_groups   : (int16) for each group (1st row), the type (2nd row) (-1 artifact, 0 unclustered, 1 multi-unit, 2 single unit)
%%%     -types_orig           :        2 x n_groups   : (int16) same as above, before manual clustering (not used)
%%%     -artifacts            :        2 x n_classes  : (int64) for each class (1st row), whether it is artifact (1) or not (0) (2nd row)
%%%     -artifacts_prematch   :        2 x n_classes  : (int64) same as above, before manual clustering (not used)
%%%     -groups               :        2 x n_classes  : (int16) for each class (1st row), which group (2nd row)
%%%     -groups_orig          :        2 x n_classes  : (int16) same as above, before manual clustering (not used)
%%%     -classes              : (n_spikes-skipped) x 1: (uint16) class number for each spike
%%%     -distance             : (n_spikes-skipped) x 1: (single) template matching distances
%%%     -index                : (n_spikes-skipped) x 1: (uint32) python index (0:length(times)-1) of the spike times corresponding to spike (some skipped)
%%%     -matches              : (n_spikes-skipped) x 1: (int8) 0 (SPC), 1 (1st template matching), 2 (second template matching)

%%% List the folders with combinato results for each channel
channel_folders = dir(combinato_directory);
channel_folders = channel_folders([channel_folders.isdir]);
channel_folders = {channel_folders.name}';
channel_folders = channel_folders(contains(channel_folders, 'NS6'));

%%% List polarity signs
polarity_signs = [{'neg'}; {'pos'}];

%%% Loop through combinations of channel, options, polarity sign
switch clustering_mode
    case 'standard'
        parameter_options = 0; %%% Default options
    case 'grid-search'
        parameter_options = (0:16)'; %%% Default options plus grid parameters
end

[Ax, Bx, Cx] = ndgrid(1:numel(polarity_signs), 1:numel(parameter_options), 1:numel(channel_folders));
T = table;
T.channel_folder = channel_folders(Cx(:));
T.parameter_option = parameter_options(Bx(:));
T.polarity_sign = polarity_signs(Ax(:));

channel_folders = T.channel_folder;
parameter_options = [T.parameter_option];
polarity_signs = T.polarity_sign;
n_combos = height(T); clear T

%%% Empty arrays to hold on to structures
channels = [];
neurons = [];
classes = [];

for idx = 1:n_combos
    
    this_folder = channel_folders{idx};
    polarity_sign = polarity_signs{idx};
    parameter_option = parameter_options(idx);
       
    channel_number = regexp(this_folder, '[0-9][0-9][0-9]', 'match');
    noise_standard_deviation = noise_info(contains(noise_info.file_name, channel_number{1}), :).noise_standard_deviation;
    channel_number = str2double(channel_number{1});
    
    spike_data_file = dir(fullfile(combinato_directory, this_folder, 'data_*.h5'));
    sorting_results_file = dir(fullfile(combinato_directory, this_folder, sprintf('sort_%s_options%02d', polarity_sign, parameter_option), 'sort_cat.h5'));
    
    if length(sorting_results_file) == 1 && length(spike_data_file) == 1
        
        %%% Read variables from spike data file and sorting results file
        spike_data_file = fullfile(spike_data_file.folder, spike_data_file.name);
        sorting_results_file = fullfile(sorting_results_file.folder, sorting_results_file.name);
    
        spikes = struct('waveforms', [], 'times', [], 'artifacts', []);
        spikes.waveforms = h5read(spike_data_file, sprintf('/%s/spikes', polarity_sign));
        spikes.times = h5read(spike_data_file, sprintf('/%s/times', polarity_sign));
        spikes.artifacts = h5read(spike_data_file, sprintf('/%s/artifacts', polarity_sign));
        
        sorting_results = struct('types', [], 'artifacts', [], 'groups', [], 'classes', [], 'index', [], 'matches', [], 'distances', []);
        sorting_results.types = h5read(sorting_results_file, '/types');
        sorting_results.artifacts = h5read(sorting_results_file, '/artifacts');
        sorting_results.groups = h5read(sorting_results_file, '/groups');
        sorting_results.classes = h5read(sorting_results_file, '/classes');
        sorting_results.index = h5read(sorting_results_file, '/index');
        sorting_results.matches = h5read(sorting_results_file, '/matches');
        sorting_results.distances = h5read(sorting_results_file, '/distance');
        
        
        %%% Make channel_info structure for use in building channels, neurons and classes structures
        channel_info = make_channel_structure(subject, folder, bank, channel_number, parameter_option, polarity_sign, file_length, sampling_rate, noise_standard_deviation);
        
        %%% Channel info logged, artifacts removed, spikes ordered by group then index, with group number (1) start (2) and end (3) indices provided.
        [this_channel, times, ISIs, waveforms, matches, distances, spike_classes, class_groups, group_indices] = transform_spike_data(channel_info, spikes, sorting_results);
        clear spikes sorting_results
        
        %%% Development: Could save artifact free and group-ordered variables in a lighter format, 
        %%% but for now will keep data in combinato format and only log information of each hierarchical
        %%% level for use in accessing and refining data in analysis.
        
        if this_channel.n_classes > 0 && this_channel.n_groups > 0
            
            %%% Classes structure created, class info logged, class silhouette
            %%% scores calculated. Classes structure created before neuron
            %%% structure to pass as an argument silhouette scores
            
            these_classes = make_class_structure(channel_info, this_channel, ISIs, waveforms, matches, distances, spike_classes, class_groups);
            
            %%% Neurons (group) structure created, group info logged
            [this_channel, these_neurons] = make_neuron_structure(channel_info, this_channel, ISIs, waveforms, matches, distances, spike_classes, class_groups, group_indices);
            
            %%% Concatenating structures associated to a session's bank
            channels = [channels; this_channel];
            neurons = [neurons; these_neurons];
            classes = [classes; these_classes];
            
            clear this_channel these_neurons these_classes
        end
        
        clear times ISIs waveforms matches distances spike_classes group_indices classes_mean_waveforms

    end
    
end

save(fullfile(save_directory, 'channels.mat'), 'channels');
save(fullfile(save_directory, 'neurons.mat'), 'neurons');
save(fullfile(save_directory, 'classes.mat'), 'classes');

end


function channel_info = make_channel_structure(subject, folder, bank, channel_number, parameter_option, polarity_sign, file_length, sampling_rate, noise_standard_deviation)
channel_info = struct;
channel_info.subject = subject;
channel_info.folder = folder;
channel_info.bank = bank;
channel_info.channel_number = channel_number;
channel_info.parameter_option = parameter_option;
channel_info.is_neg = strcmp(polarity_sign, 'neg');
channel_info.file_length = file_length;
channel_info.sampling_rate = sampling_rate;
channel_info.noise_standard_deviation = noise_standard_deviation;
end


function [channel, times, ISIs, waveforms, spike_matches, spike_distances, spike_classes, class_groups, group_indices] = transform_spike_data(channel_info, spikes, sorting_results)

old_times = spikes.times; 
n_detections = length(old_times);
old_waveforms = double(spikes.waveforms);
spike_artifacts = spikes.artifacts;
spike_artifacts = find(spike_artifacts);
group_types = double(sorting_results.types);
class_groups = double(sorting_results.groups);
class_artifacts = double(sorting_results.artifacts);
class_artifacts = class_artifacts(1, class_artifacts(2, :)==1);
spike_classes = double(sorting_results.classes);
spike_indices = double(sorting_results.index);
spike_matches = double(sorting_results.matches);
spike_distances = double(sorting_results.distances);
n_skipped = size(old_waveforms, 2) - length(spike_indices);
bad_types = group_types(2, :)<1;
bad_groups = group_types(1, bad_types);
bad_classes = class_groups(1, ismember(class_groups(2, :), bad_groups));
bad_classes = unique([bad_classes, class_artifacts]);
bad_indices = spike_indices(ismember(spike_classes, bad_classes));
bad_indices = unique([bad_indices;(spike_artifacts-1)]); 

%%% Combinato uses python indices, so subtract 1 from artifacts to find good indices

%%% Doing it this way because some spike indices are skipped
good_spikes = ~ismember(spike_indices, bad_indices);
good_indices = spike_indices(good_spikes);
spike_classes = spike_classes(good_spikes);
spike_matches = spike_matches(good_spikes);
spike_distances = spike_distances(good_spikes);
class_groups = class_groups(:, ismember(class_groups(1, :), unique(spike_classes)));
table_groups = arrayfun(@(x) class_groups(2, class_groups(1, :)==x), spike_classes);

%%% Make table to sort spike data by group then index
index_table = table;
index_table.index = good_indices;
index_table.class = spike_classes;
index_table.group = table_groups;
index_table.distance = spike_distances;
index_table.match = spike_matches;
index_table = sortrows(index_table, {'group', 'index'}, 'ascend');

spike_classes = index_table.class;
spike_distances = index_table.distance;
spike_matches = index_table.match;
table_groups = index_table.group; 
good_indices = index_table.index;
clear index_table;

good_indices = good_indices + 1; %%% Add + 1 to keep only the good spikes.

times = old_times(good_indices); clear old_times %%% Sorted by group and index
waveforms = old_waveforms(:, good_indices);
waveforms = waveforms';
ISIs = diff(times); %%% Interspike intervals
ISIs(ISIs < 0)= NaN;

good_classes = unique(spike_classes, 'stable');
good_groups = unique(table_groups, 'stable');

group_indices = arrayfun(@(x) [x, find(table_groups == x, 1, 'first'), find(table_groups == x, 1, 'last')], good_groups, 'UniformOutput', false);
group_indices = cell2mat(group_indices);

%%% Create channel structure with info about clustering results
channel = channel_info;
channel.n_detections = n_detections;
channel.n_groups = length(good_groups);
channel.n_bad_groups = length(bad_groups);
channel.n_classes = length(good_classes);
channel.n_bad_classes = length(bad_classes);
channel.n_spikes = length(times);
channel.n_bad_spikes = length(bad_indices);
channel.n_skipped = n_skipped;
channel.n_features = size(waveforms, 2);
channel.n_spc = sum(spike_matches == 0);
channel.n_1st = sum(spike_matches == 1);
channel.n_2nd = sum(spike_matches == 2);

end


function classes = make_class_structure(channel_info, channel, ISIs, waveforms, matches, distances, spike_classes, class_groups)
%%% Get channel info
n_classes = channel.n_classes;
n_groups = channel.n_groups;
n_features = channel.n_features;
n_detections = channel.n_detections;
is_neg = channel.is_neg;

%%% Define line noise interspike intervals in ms, floored for indexing in histogram
harmonics_ISIs = floor(1000./(60.*(1:5))); %%%in ms the expected ISIs of 60Hz to 5th harmonic
pseudo_subharmonics_ISIs = floor(1000./(60./(2:15))); %%%in ms ISIs of seeming subharmonics occuring in multiples of 16.66ms
line_noise_ISIs = [harmonics_ISIs, pseudo_subharmonics_ISIs];

%%% Create array with individual spike group labels
unique_classes = class_groups(1, :)';

%%% Start classes structure from channel_info
classes = repmat(channel_info, n_classes, 1);

%%% Loop through classes
for kdx = 1:n_classes
    
    %%% Get class and group number
    this_class = unique_classes(kdx);
    this_group = class_groups(2, unique_classes==this_class);
    
    %%% Find indices of spikes in the class
    spike_indices = find(spike_classes == this_class);
    n_spikes = length(spike_indices); %%%n_spikes in this class
    
    %%% Get spike data
    class_matches = matches(spike_indices);
    class_distances = distances(spike_indices);
    
    %%% For ISI indices subtract one and remove 0 index (first spike not counted)
    ISI_indices = spike_indices - 1;
    ISI_indices = ISI_indices(ISI_indices ~= 0);
    class_ISIs = ISIs(ISI_indices); %%% ISIs are as calculated within entire group.
 
    %%% Calculate percentage of ISIs below 3 ms and contributing to line noise
    p_sub3ms = mean(class_ISIs < 3); %%% Proportion of spikes in this class that contribute to < 3 ms in group
    histogram_ISIs = histcounts(class_ISIs, 1:256); %%% Not counting 0 to 1 for indices to match edges
    p_60Hz = sum(histogram_ISIs(line_noise_ISIs))/n_spikes; %%% Proportion of spikes in this class that contribute to line noise ISIs
    
    if n_spikes > n_features
        
        %%% Calculate the Mahalanobis distances of spikes in class
        centroid = mean(waveforms(spike_indices, :), 1);
        inverse_covariance_matrix = inv(cov(waveforms(spike_indices, :)));
        differences = waveforms(spike_indices, :) - centroid;
        intra_mahalanobis = sqrt(sum(differences * inverse_covariance_matrix .* differences, 2));
        clear differences;
        
        %%% Class mean dispersion is mean of squared distances
        m_dispersion = mean(intra_mahalanobis.^2);
        
        %%% Gathering the mean and std Mahalanobis as measures of cohesion
        m_mahalanobis = mean(intra_mahalanobis);
        std_mahalanobis = std(intra_mahalanobis);
        clear intra_mahalanobis
        
        %%% If there is more than one class get isolation distances with respects
        %%% to spikes in other classes of same group (cohesion) or other groups (separation).
        if n_classes > 1
            
            %%% Define other classes
            other_classes = unique_classes(unique_classes ~= this_class);
            
            %%% Define other groups per other class
            other_class_groups = class_groups(2, unique_classes~=this_class);
            
            %%% In the same group
            intra_group_classes = other_classes(other_class_groups==this_group);
            
            if ~isempty(intra_group_classes)
                
                %%% Get isolation distance of all other spikes in group
                other_class_spikes = find(ismember(spike_classes, intra_group_classes));
                
                %%% If there are more spikes in class than outside class will use
                %%% the maximum distance determine by number of other spikes
                n_min = min([length(spike_indices), length(other_class_spikes)]);
                
                %%% Calculate Mahalanobis distances of other spikes
                differences = waveforms(other_class_spikes, :) - centroid;
                extra_mahalanobis = sqrt(sum(differences * inverse_covariance_matrix .* differences, 2));
                clear differences
                
                %%% Sort them ascending and select N_min'th distance as isoD
                extra_mahalanobis = sort(extra_mahalanobis, 'ascend');
                intragroup_isolation_distance = extra_mahalanobis(n_min);
                clear extra_mahalanobis
            else
                intragroup_isolation_distance = NaN;
            end
            
            %%% If there are other groups
            if n_groups > 1
            
                %%% Measure isolation distance based on spikes from other groups
                extra_group_classes = other_classes(other_class_groups~=this_group);
                
                other_class_spikes = find(ismember(spike_classes, extra_group_classes));
                n_min = min([length(spike_indices), length(other_class_spikes)]);
                
                differences = waveforms(other_class_spikes, :) - centroid;
                extra_mahalanobis = sqrt(sum(differences * inverse_covariance_matrix .* differences, 2));
                clear differences
                
                extra_mahalanobis = sort(extra_mahalanobis, 'ascend');
                extragroup_isolation_distance = extra_mahalanobis(n_min);
                clear extra_mahalanobis
            else
                extragroup_isolation_distance = NaN;
            end
            
            %%% For L-ratio calculate all squared Mahalanobis distances in set
            %%% with respect to class's covariance matrix
            differences = waveforms - centroid;
            all_mahalanobis_squared = sum(differences * inverse_covariance_matrix .* differences, 2);
            
            %%% 1-chicdf is probability each spike should belong to class
            P = 1 - chi2cdf(all_mahalanobis_squared, n_features);
            clear all_mahalanobis_squared
            
            %%% Get the probabilities for spikes outside the class, assuming all
            %%% spikes in class have the greatest probabilities
            P_extra = P(length(spike_indices) + 1:end);
            
            %%% L-ratio is sum of extraclass probabilites divided by class length
            L_ratio = sum(P_extra)/length(spike_indices);
            
        else
            %%% If not more than one class these metrics are NaN
            intragroup_isolation_distance = NaN;
            extragroup_isolation_distance = NaN;
            L_ratio = NaN;
        end
    else
        m_dispersion = NaN;
        m_mahalanobis = NaN;
        std_mahalanobis = NaN;
        intragroup_isolation_distance = NaN;
        extragroup_isolation_distance = NaN;
        L_ratio = NaN;
    end
    
    %%% Get class absolute waveforms for many other metrics
    class_waveforms = waveforms(spike_indices, :);
    if is_neg
        class_waveforms = class_waveforms * -1; %%% In order to get absolute peak values
    end
    
    %%% Log class info
    
    %%% Class and group numbers
    classes(kdx).class = this_class;
    classes(kdx).group = this_group;
    
    %%% The total number of spikes in the class, and as a faction from detected in the channel recording
    classes(kdx).n_spikes = n_spikes;                           
    classes(kdx).f_spikes = n_spikes/n_detections;
    
    %%% The percentage of spikes in this class that contribute to sub3ms and line noise ISIs
    classes(kdx).p_sub3ms = p_sub3ms;
    classes(kdx).p_60Hz = p_60Hz;
    
    %%% Percentages in this class that were matched through SPC, 1st template or 2nd template matching
    classes(kdx).p_spc = sum(class_matches == 0)/length(class_matches);    
    classes(kdx).p_1st = sum(class_matches == 1)/length(class_matches);
    classes(kdx).p_2nd = sum(class_matches == 2)/length(class_matches);
    
    %%% Mean and standard deviation of the class's spike amplitudes in uV (absolute value).
    classes(kdx).m_amp_uV = mean(max(class_waveforms, [], 1));              
    classes(kdx).std_amp_uV = std(max(class_waveforms, [], 1));
    
    %%% Mean and standard deviation of average voltage of class's waveforms
    classes(kdx).m_voltage = mean(mean(class_waveforms, 1));                
    classes(kdx).std_voltage = std(mean(class_waveforms, 1));               
    
    %%% Mean and standard deviation of all matching distances
    classes(kdx).m_distances = mean(class_distances);                      
    classes(kdx).std_distances = std(class_distances);
    
    %%% Calculated dispersion
    classes(kdx).m_dispersion = m_dispersion;                                  
    
    %%% Mean and standard deviation of mahalanobi distances of spikes in class to its centroid
    classes(kdx).m_mahalanobis = m_mahalanobis;                     
    classes(kdx).std_mahalanobis = std_mahalanobis;
    
    %%% Isolation distance of class to other classes in group
    classes(kdx).intragroup_isolation_distance = intragroup_isolation_distance;
    
    %%% Isolation distance of class to other classes outside group
    classes(kdx).extragroup_isolation_distance = extragroup_isolation_distance;
    
    %%% L-ratio
    classes(kdx).L_ratio = L_ratio;
end

end


function [channel, neurons] = make_neuron_structure(channel_info, channel, ISIs, waveforms, matches, distances, spike_classes, class_groups, group_indices)
%%% Get channel info
n_groups = channel.n_groups;
n_detections = channel.n_detections;
n_features = channel.n_features;
is_neg = channel.is_neg;
file_length = channel.file_length;
sampling_rate = channel.sampling_rate;
duration_seconds = file_length/sampling_rate;

%%% Define line noise interspike intervals in ms, floored for indexing in histogram
harmonics_ISIs = floor(1000./(60.*(1:5))); %%% in ms the expected ISIs of 60Hz to 5th harmonic
pseudo_subharmonics_ISIs = floor(1000./(60./(2:15))); %%% in ms ISIs of seeming subharmonics occuring in multiples of 16.66ms
line_noise_ISIs = [harmonics_ISIs, pseudo_subharmonics_ISIs];
not_noise_ISIs = ~ismember(1:255, line_noise_ISIs);

%%% Start neurons (group) structure from channel_info
neurons = repmat(channel_info, n_groups, 1);

%%% Make a spike group label array
unique_classes = class_groups(1, :);
spike_groups = arrayfun(@(x) class_groups(2, unique_classes == x), spike_classes);
unique_groups = group_indices(:, 1);

if n_groups > 1
    %%% Preallocate array to hold values for Davies-Bouldin Index
    DBI_D = zeros(n_groups, 1);
    
    %%% Define global centroid for Calinski-Harabasz Index
    global_centroid = mean(waveforms, 1);
    
    %%% Preallocate array to hold Calinski Bs
    calinski_B = zeros(n_groups, 1);
    calinski_W = zeros(n_groups, 1);
    
    n_total = channel.n_spikes;
end

%%% Loop through groups
for kdx = 1:n_groups
    
    %%% Get indices of spikes for the group
    this_group = unique_groups(kdx);
    spike_indices = find(spike_groups == this_group);
    
    %%% Get number of spikes and firing rate
    n_spikes = length(spike_indices);
    firing_rate = n_spikes/duration_seconds;
    
    %%% Determine how many classes in group
    group_classes = class_groups(2, :)==this_group;
    n_classes = sum(group_classes);
    
    %%% Get spike data
    group_spike_matches = matches(spike_indices);
    group_spike_distances = distances(spike_indices);
    
    %%% For ISIs subtract 1 from spike_indices and remove 0 index (first spike)
    ISI_idx = spike_indices - 1;
    ISI_idx = ISI_idx(ISI_idx~=0);
    group_ISIs = ISIs(ISI_idx);
    
    %%% Calculate percentage of ISIs below 3 ms and contributing to line noise
    p_sub3ms = mean(group_ISIs<3); %%%Percentage of ISIs in group < 3ms
    histogram_ISIs = histcounts(group_ISIs, 1:256); %%% Not counting 0 to 1 for indices to match edges
    p_60Hz = sum(histogram_ISIs(line_noise_ISIs))/n_spikes; %%% Percentage of ISIs in line noise ISIs
    
    %%% Calculate ratio of mean line_noise counts to not noise counts
    isi_SNR = mean(histogram_ISIs(not_noise_ISIs))/mean(histogram_ISIs(line_noise_ISIs));
    
    if n_spikes > n_features %%% Having less observations than features can lead to errors due to issues calculating the inverse of the covariance matrix
        
        %%% Calculate the Mahalanobis distances of spikes in group
        centroid = mean(waveforms(spike_indices, :), 1);
        inverse_covariance_matrix = inv(cov(waveforms(spike_indices, :)));
        differences = waveforms(spike_indices, :)-centroid;
        intra_mahalanobis = sqrt(sum(differences * inverse_covariance_matrix .* differences, 2));
        clear differences;
        
        %%% Group mean dispersion is mean of squared distances
        m_dispersion = mean(intra_mahalanobis.^2);
        
        %%% Gathering the mean and std Mahalanobis as measures of cohesion
        m_mahalanobis = mean(intra_mahalanobis);
        std_mahalanobis = std(intra_mahalanobis);
        
        %%% If there is more than one group get quality metrics that are based on intercluster (separation) metrics.
        if n_groups > 1
            
            %%%Define other groups
            other_groups = unique_groups(unique_groups~=this_group);
            
            %%% Get isolation distance of all other spikes outside group
            other_group_spikes = find(ismember(spike_groups, other_groups));
            
            %%% If there are more spikes in group than outside group will use the maximum distance determined by number of other spikes
            n_min = min([length(spike_indices), length(other_group_spikes)]);
            
            %%% Calculate Mahalanobis distances of other spikes
            differences = waveforms(other_group_spikes, :)-centroid;
            extra_mahalanobis = sqrt(sum(differences * inverse_covariance_matrix .* differences, 2));
            clear differences
            
            %%% Sort them ascending and select N_min'th distance as isolation_distance
            extra_mahalanobis = sort(extra_mahalanobis, 'ascend');
            group_isolation_distance = extra_mahalanobis(n_min);
            clear extra_mahalanobis
            
            %%% For L-ratio calculate all squared Mahalanobis distances in set with respect to group's covariance matrix
            differences = waveforms - centroid;
            all_mahalanobis_squared = sum(differences * inverse_covariance_matrix .* differences, 2);
            
            %%% 1-chicdf is probability each spike should belong to group
            P = 1 - chi2cdf(all_mahalanobis_squared, n_features);
            clear all_mahalanobis_squared
            
            %%% Get the probabilities for spikes outside the group, assuming all spikes in group have the greatest probabilities
            P_extra = P(length(spike_indices)+1:end);
            
            %%% L-ratio is sum of extragroup probabilites divided by group length
            L_ratio = sum(P_extra)/length(spike_indices);
            
            %%% For Dunn's index, calculating as mean Mahalanobis distance of
            %%% cluster divided by the minimum intergroup distance (Mahalanobis
            %%% distances of other group centroids with resepect to group's
            %%% covariance matrix
            other_group_spikes = arrayfun(@(x) find(spike_groups == x), other_groups, 'UniformOutput', false);
            other_centroids = cell2mat(cellfun(@(x) mean(waveforms(x, :), 1), other_group_spikes, 'UniformOutput', false));
            differences = other_centroids - centroid;
            other_group_mahalanobiss = sqrt(sum(differences * inverse_covariance_matrix .* differences, 2));
            dunn_index = m_mahalanobis/min(other_group_mahalanobiss);
            
            %%% For Davies-Bouldin Index (at channel level) calculating the
            %%% mean Mahal distances of other spikes by group, then R(i,j) would be
            %%% equal to (m_mahalanobisi + m_mahalanobisj)/other_group_mahalanobiss(j)
            other_m_mahalanobis = cell2mat(cellfun(@(x) mean(sqrt(sum((waveforms(x, :) - centroid) * inverse_covariance_matrix .* (waveforms(x, :) - centroid), 2))), other_group_spikes, 'UniformOutput', false));
            DBI_D(kdx) = max(arrayfun(@(x) (m_mahalanobis + other_m_mahalanobis(x)) / other_group_mahalanobiss(x), 1:length(other_groups)));
            
            %%% For Calinski-Harabasz Index using B: sum of squared Mahalanobis distances
            %%% between group centroids and the global centroid (times the number
            %%% of spikes in group), and W the sum of squared Mahalanobis
            %%% distances for spikes in the group to their centroid
            calinski_W(kdx) = n_spikes*m_dispersion;
            difference = centroid-global_centroid;
            calinski_B(kdx) = n_spikes*sum(difference*inverse_covariance_matrix.*difference, 2);
            
            %%% For Silhouette scores a(i) will be Mahalanobis distance of spike i
            %%% in this group to centroid with respect to group covariance matrix
            %%% and b(i) will be min Mahalanobis distance of spike i to other group
            %%% centroid with respect to other groups covariance matrices
            other_cov_mats = cellfun(@(x) inv(cov(waveforms(x, :))), other_group_spikes, 'UniformOutput', false);
            silhouette_a = intra_mahalanobis; clear intra_mahalanobis
            silhouette_b = arrayfun(@(x) sqrt(sum((waveforms(spike_indices, :) - other_centroids(x, :)) * other_cov_mats{x} .* (waveforms(spike_indices, :) - other_centroids(x, :)), 2)), 1:length(other_groups), 'UniformOutput', false);
            silhouette_b = cell2mat(silhouette_b);
            silhouette_b = min(silhouette_b, [], 2);
            silhouette_scores = (silhouette_a-silhouette_b)./max([silhouette_a, silhouette_b], [], 2);
            silhouette_score = mean(silhouette_scores);
        else
            %%% If not more than one group, these metrics are NaN
            group_isolation_distance = NaN;
            L_ratio = NaN;
            dunn_index = NaN;
            silhouette_score = NaN;
        end
    else
        m_dispersion = NaN;
        m_mahalanobis = NaN;
        std_mahalanobis = NaN;
        group_isolation_distance = NaN;
        L_ratio = NaN;
        dunn_index = NaN;
        silhouette_score = NaN;
    end
    
    group_spikes = waveforms(spike_indices, :);
    if is_neg
        group_spikes = group_spikes * -1; %%% In order to get absolute peak values
    end
    centroid = mean(group_spikes, 1);
    group_std_spikes = std(group_spikes, [], 1);
    
    %%% Get mean waveform peak, width, peak std and avg std of last 25% of samples
    mean_peak_uV = max(centroid);
    mean_peak_idx = find(centroid == mean_peak_uV, 1, 'first');
    mean_trough = min(centroid(1:mean_peak_idx-1));
    mean_trough_idx = find(centroid(1:mean_peak_idx-1) == mean_trough, 1, 'first');
    peak_width =  mean_peak_idx - mean_trough_idx + 1;
    std_peak = std(group_spikes(:, mean_peak_idx));
    std_last25 = mean(group_std_spikes(round(n_features*0.75):end));
    inflection_points = find(diff(diff(centroid))==0);
    inflection_points = inflection_points(inflection_points > mean_peak_idx - 2);
    if ~isempty(inflection_points)
        inflection_values = [mean_peak_uV, centroid(inflection_points)];
        inflection_magnitudes = abs(diff(inflection_values));
        num_inflection_violations = sum(inflection_magnitudes > (0.2*abs(mean_peak_uV)));
    else
        num_inflection_violations = 0;
    end
    
    %%% Get group grade based on criteria
    grade = grade_group(firing_rate, p_sub3ms, p_60Hz, isi_SNR, mean_peak_uV, std_last25, std_peak, num_inflection_violations);
    
    %%% Log class info
    
    %%% Group number
    neurons(kdx).group = this_group;
    
    %%% Grade divided in 3 bins to filter out a certain type
    neurons(kdx).is_SUA = strcmp(grade, 'SUA');
    neurons(kdx).is_MUA = strcmp(grade, 'MUA');
    neurons(kdx).is_noise = strcmp(grade, 'noise');
    
    %%% The total number of spikes in the class, and as a faction from detected in the channel recording
    neurons(kdx).n_spikes = n_spikes;                           
    neurons(kdx).f_spikes = n_spikes/n_detections;
    neurons(kdx).n_classes = n_classes;
    
    %%% The mean firing rate in Hz
    neurons(kdx).mean_FR = firing_rate;
    neurons(kdx).mean_peak_uV = mean_peak_uV;
    neurons(kdx).peak_width = peak_width;
    
    %%% The percentage of spikes in this class that contribute to sub3ms and line noise ISIs
    neurons(kdx).p_sub3ms = p_sub3ms;
    neurons(kdx).p_60Hz = p_60Hz;
    neurons(kdx).isi_SNR = isi_SNR;
    
    %%% Percentages in this class that were matched through SPC, 1st template or 2nd template matching
    neurons(kdx).p_spc = sum(group_spike_matches == 0)/length(group_spike_matches);    
    neurons(kdx).p_1st = sum(group_spike_matches == 1)/length(group_spike_matches);
    neurons(kdx).p_2nd = sum(group_spike_matches == 2)/length(group_spike_matches);
    
    %%% Mean and standard deviation of the group's spike amplitudes in uV (absolute value).
    neurons(kdx).m_amp_uV = mean(max(group_spikes, [], 1));              
    neurons(kdx).std_amp_uV = std(max(group_spikes, [], 1));
    
    %%% Mean and standard deviation of average voltage of class's waveforms
    neurons(kdx).m_voltage = mean(mean(group_spikes, 1));                
    neurons(kdx).std_voltage = std(mean(group_spikes, 1));               
    
    %%% Mean and standard deviation of all matching distances
    neurons(kdx).m_distances = mean(group_spike_distances);                      
    neurons(kdx).std_distances = std(group_spike_distances);
    
    %%% Calculated dispersion
    neurons(kdx).dispersion = m_dispersion;                                  
    
    %%% Mean and standard deviation of mahalanobi distances of spikes in group to its centroid
    neurons(kdx).m_mahalanobis = m_mahalanobis;                     
    neurons(kdx).std_mahalanobis = std_mahalanobis;
    
    %%% Isolation distance of group to other groups in channel
    neurons(kdx).isolation_distance = group_isolation_distance;
    
    %%% L_ratio
    neurons(kdx).L_ratio = L_ratio;
    
    %%% Dunn Index
    neurons(kdx).dunn_index = dunn_index;
    
    %%% Average silhouette scores of classes in group
    neurons(kdx).silhouette = silhouette_score;
end

%%% Get channel totals
channel.n_sua = sum([neurons.is_SUA]);
channel.n_mua = sum([neurons.is_MUA]);
channel.n_noise = sum([neurons.is_noise]);
channel.isi_SNR = mean([neurons.isi_SNR]);

%%% Get channel quality metrics encompassing group measures
if n_groups > 1 && n_spikes > n_features
    channel.DBI = mean(DBI_D);
    channel.CHI = (sum(calinski_B)/sum(calinski_W))*((n_total-n_groups)/(n_groups-1));
    channel.m_silhouette = mean([neurons.silhouette]);
else
    channel.DBI = NaN;
    channel.CHI = NaN;
    channel.m_silhouette = NaN;
end
    
end


function grade = grade_group(firing_rate, p_sub3ms, p_60Hz, isi_SNR, mean_peak_uV, std_last25, std_peak, num_inflection_violations)
%%% Grading criteria adapted from doi:10.1088/1741-2560/10/1/016001
%%% Substituted PSD criteria for other line noise features

grade = 'noise'; %%% Until proven otherwise

%%% Criterion 1: Firing rate above a threshold
fr_sua = firing_rate > 0.05;
fr_mua = firing_rate > 0.05;

%%% Criterion 2: Refractory period violation below a threshold
sub3ms_sua = p_sub3ms < 0.05;
sub3ms_mua = p_sub3ms < 0.1;

%%% Criterion 3: Ratio of peak amplitude to standard deviation of last 25% of samples
peak_ratio = mean_peak_uV/std_last25;
peak_ratio_sua = peak_ratio > 1;
peak_ratio_mua = peak_ratio > 0.5;

%%% Criterion 4: Ratio of peak std to stds of last 25% of samples
std_ratio = std_peak/std_last25;
std_ratio_sua = std_ratio > 0.33;
std_ratio_mua = std_ratio > 0.33;

%%% Criterion 5: Number of inflection magnitude violations
inflections_sua = num_inflection_violations < 3;
inflections_mua = num_inflection_violations <= 3;

%%% Must meet all criteria to be graded mua or sua
if fr_mua && sub3ms_mua && peak_ratio_mua && std_ratio_mua && inflections_mua
    grade = 'MUA';
end
if fr_sua && sub3ms_sua && peak_ratio_sua && std_ratio_sua && inflections_sua
    grade = 'SUA';
end

%%% Adding this to criteria. If exceeds certain noise levels then grade as noise instead
if p_60Hz > 0.2 || isi_SNR < 0.8
    grade = 'noise';
end

end
