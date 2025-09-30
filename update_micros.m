%%% Update Micros Function
%%%
%%% This is the master function for the micros pipeline, outlining the steps
%%% carried out in processing recording files in order to extract and cluster
%%% neural spiking data with the objective of acquiring time series of
%%% single unit activity.
%%%
%%% Function will go through 
%%%
%%% Nodes with at least 256GB of RAM should be used. 32GB of RAM will be insufficient
%%% for pipeline step #3 denoising using 'zapline' method without spliting recording in chunks.
%%% For all other steps, memory usage has been measured to select number of workers
%%% that would execute function in parallel to maximize RAM usage.
%%%
%%% Function can be executed from bash script using the SLURM to process multiple 
%%% recordings simultaneously.
%%% The first part of the code serves to organize execution of function in multiple nodes 
%%% so that a single one of them will update the table tracking progress for each recording.
%%%
%%% Function can also be executed from editor using desired inputs.

function update_micros(varargin)
if isempty(varargin)                                                 
    
    %%% To run from editor edit these values
    
    using_SLURM = false;                                             %%% Whether running pipeline from multiple nodes. Always false running from editor.
    set_threshold = 0;                                               %%% Threshold time in seconds for an individual node to start. 0 running from editor.
    ellapsed_time = 0;                                               %%% Any process ID files created within an amount of time lesser than this will be deleted.
    start_time = datestr(datetime('now'), 'yyyymmddHHMMSS');         %%% Time of function execution
    root_directory = '/path/to/micros_pipeline/parent_directory';    %%% Parent directory containing 'micros_pipeline' and 'micros_database' folders.
    referencing_method = 'reref_func';                               %%% Options: 'reref_func'
    referencing_mode = 'signal_power';                               %%% Options: 'common_average', 'signal_power', 'none'
    filtering_method = 'zapline';                                    %%% Options: 'zapline'
    clustering_method = 'combinato_func';                            %%% Options: 'combinato_func'
    clustering_mode = 'standard';                                    %%% Options: 'standard', 'grid-search'
    
else

    %%% Otherwise this is the order of function inputs bash script should follow
    
    using_SLURM = true;
    set_threshold = 300;
    ellapsed_time = -15;
    start_time = varargin{1};
    root_directory = varargin{2};
    referencing_method = varargin{3};
    referencing_mode = varargin{4};
    filtering_method = varargin{5};
    clustering_method = varargin{6};
    clustering_mode = varargin{7};
    
end

warning off

%%% For each step in the pipeline, selected number of workers for parallelized processing to maximize 256GB RAM usage
n_workers_split = 8;
n_workers_reref = 3;
n_workers_rescale = 8;
n_workers_5to8 = 10;
n_workers_lfp = 8;
n_workers_modal = 5;
n_workers_align = 18;

%%% Directory for process IDs
process_ID_directory = fullfile(root_directory, 'micros_pipeline/process_files/process_IDs');

time_threshold = set_threshold; %%% In seconds, the permitted ellapsed time
                                %%% since start_time for a node to begin
                                %%% the pipeline steps. If a node exceeds this time,
                                %%% it will quit the code to avoid errors. In the
                                %%% first step, the last node to finish will pause
                                %%% for this amount of time to allow all other nodes
                                %%% that started within this interval to come into
                                %%% sync.
                                
start_datenum = datenum(start_time, 'yyyymmddHHMMSS'); %%% In days with fraction

%%% Pipeline initialization steps:

%%% 1) First delete all process IDs that do not have the same start_time. These
%%% are the process IDs from a previous SLURM attempt that failed.
process_ID_files = dir(fullfile(process_ID_directory, '*_process_ID_*'));
if ~isempty(process_ID_files)
    file_datenums = [process_ID_files.datenum];
    time_difference = (file_datenums - start_datenum) * 24 * 60; % In minutes    
    bad_files = process_ID_files(time_difference < ellapsed_time); %15 minutes too long and reasonable to start new job?
    if ~isempty(bad_files)
        bad_files = fullfile({bad_files.folder},{bad_files.name});
        for idx = 1:length(bad_files)
            bad_file = bad_files{idx};
            delete(bad_file);
        end
    end
end

%%% 2) Create a new file with the process ID for this node labeled with start_time
process_ID = feature('getpid');
process_ID_file = fullfile(process_ID_directory, sprintf('%s_process_ID_%d', start_time, process_ID));
file_ID = fopen(process_ID_file, 'w'); 
fclose(file_ID);
fprintf('%s - Node %d starting.\n', datestr(now), process_ID);

%%% 3) Small pause for process ID files to be written before continuing
pause(1)

%%% 4) Check directory to get all process IDs of files written so far
process_ID_files = dir(fullfile(process_ID_directory, '*_process_ID_*'));
process_ID_datenums = [process_ID_files.datenum];
process_ID_files = {process_ID_files.name};  
process_ID_files = strrep(process_ID_numbers, start_time, '');
process_ID_numbers = strrep(process_ID_numbers, '_process_ID_', '');
process_ID_numbers = cellfun(@(x) str2double(x), process_ID_numbers);

%%% 5) Make node quit if time difference from start time to process ID file creation is
%%% greater than time threshold. To avoid nodes becoming stuck.
process_ID_datenum = process_ID_datenums(process_ID_numbers == process_ID);
time_difference = (process_ID_datenum - start_datenum) * 24 * 60 * 60; %%% In seconds

if using_SLURM
    if time_difference >= time_threshold
        fprintf('%s - Node %d quit because it was late.\n', datestr(now), process_ID);
        quit
    end
end

%%% 6) Select the earliest node to be the leader node. Millisecond precision.
%%% Millisecond precision could be insufficient, but only one will be chosen.
[~, sorting_indices] = sort(process_ID_datenums, 'ascend');
process_ID_numbers = process_ID_numbers(sorting_indices);
node_position = find(process_ID_numbers == process_ID);
is_leader_node = node_position == 1;

%%% 7) Separate nodes that started simultaneously by introducing a pause
%%% relative to their index.
seconds_to_pause = find(process_ID_numbers == process_ID);
pause(seconds_to_pause);

%%% 8) Set random number generator seed with process_ID for unique random numbers
rng(process_ID);

%%% Pipeline Steps:

try
    %%% Step 1. Split .ns6 recording files into separate files for each channel.
    if is_leader_node
        fprintf('%s - Leader node %d entering step 1: splitting .ns6 recording files.\n', datestr(now), process_ID);
    else
        fprintf('%s - Follower node %d entering step 1: splitting .ns6 recording files.\n', datestr(now), process_ID);
    end
    n01_split_recording_file(time_threshold, is_leader_node, n_workers_split, root_directory);
    %%% Check completion, log information into table.
    n01_split_recording_file_check(is_leader_node, root_directory);
    
    
    %%% Step 2. Rereference data: mode 'none' does highpass of 0.5Hz and removes offset,
    %%% 'common_average' does common average referencing, 'signal_power' selects the
    %%% channel with the least amount of signal power as reference.
    fprintf('%s - Node %d entering step 2: rereferencing with method: ''%s'' mode: ''%s''.\n', datestr(now), process_ID, rereferencing_method, rereferencing_mode);
    n02_rereference_signal(is_leader_node, n_workers_reref, root_directory, referencing_method, referencing_mode);
    %%% Check completion, update table.
    n02_rereference_signal_check(is_leader_node, root_directory);
    
    
    %%% Step 3. Denoise data. Methods available: 'zapline' TODO: add more methods
    fprintf('%s - Node %d entering step 3: denoising with method: ''%s''.\n',datestr(now),pid,denoising_method);
    n03_filter_noise(is_leader_node, root_directory, filtering_method);
    %%% Check completion, update table.
    n03_filter_noise_check(is_leader_node, root_directory);
    
    
    %%% Step 4. Rescale data. Combinato assumes data is in uV. 
    %%% Blackrock files need to be rescaled by a factor of 0.25 to be in uV.
    fprintf('%s - Node %d entering step 4: rescaling data to uV in int16 format.\n',datestr(now),pid);
    n04_rescale_signal(is_leader_node, n_workers_rescale, root_directory);
    %%% Check completion, update table.
    n04_rescale_signal_check(is_leader_node, root_directory);
    
    
    %%% Step 5. Convert .mat files to .hdf5 in the format accepted by Combinato.
    fprintf('%s - Node %d entering step 5: formatting data to .hdf5 supported by combinato.\n',datestr(now),pid);
    n05_convert_to_hdf5(is_leader_node, n_workers_5to8, root_directory);
    %%% Check completion, update table.
    n05_convert_to_hdf5_check(is_leader_node, root_directory);
    
    
    %%% Step 6. Run spike extraction and clustering.
    %%% Methods available: 'combinato_func'
    %%% Combinato modes available: 'standard', 'grid-search'
    fprintf('%s - Node %d entering step 6: detecting and clustering spikes with method: ''%s'' mode: ''%s''.\n',datestr(now),pid,clustering_method,clustering_mode);
    n06_cluster_signal_spikes(is_leader_node, root_directory, clustering_method, clustering_mode);
    %%% Check completion, update table.
    n06_cluster_signal_spikes_check(is_leader_node, root_directory);
    
    
    %%% Step 7. Check updates on combinato clustering files, make lists for manual review.
    fprintf('%s - Node %d entering step 7: making lists of clustering results for review.\n',datestr(now),pid);
    n07_manual_review_check(is_leader_node, root_directory);
    
    
    %%% Step 8. Collect spike and clustering data, and quality metrics at the micros bank (depth eletrode),
    %%% channel, cluster group, and cluster class levels
    fprintf('%s - Node %d entering step 8: acquiring spike clustering data and quality metrics from cluster mode: ''%s''.\n',datestr(now),pid,clustering_mode);
    n08_gather_spike_data(is_leader_node, root_directory, clustering_mode);
    %Check completion, update table.
    n08_gather_spike_data_check(is_leader_node, root_directory);
    
    
    %%% Step 9. Resample 30kHz rerefenced and denoised data to 2kHz (LFP)
    fprintf('%s - Node %d entering step 9: resampling unrescaled denoised data and saving LFPs.\n',datestr(now),pid);
    n09_get_local_field_potential(is_leader_node, n_workers_lfp, root_directory);
    %%% Check completion, update table.
    n09_get_local_field_potential_check(is_leader_node, root_directory);
    
    
    %%% Step 10. Run multiple oscillation detection algorithm (MODAL) on LFP data.
    fprintf('%s - Node %d entering step 10: running MODAL and saving results.\n',datestr(now),pid);
    n10_oscillation_detection(is_leader_node, n_workers_modal, root_directory);
    %%% Check completion, update table.
    n10_oscillation_detection_check(is_leader_node, root_directory);
    
    
    %%% Step 11. For recordings that have events files and sync pulses, align sync pulse times to event times and update events files.
    fprintf('%s - Node %d entering step 11: aligning behavioral task computer sync pulse times to Blackrock recording pulse times.\n',datestr(now),pid);
    n11_align_behavioral_data(is_leader_node, n_workers_align, root_directory);
    %%% Check completion, update table.
    n11_align_behavioral_data_check(is_leader_node, root_directory);
    
    
    %%% Post processing steps
    %%% Delete the process IDs.
    if is_leader_node
        delete(fullfile(process_ID_directory, '*'));
    end
    
catch iteration_error
    pause(time_threshold)
    fprintf('%s\n', getReport(iteration_error, 'extended', 'hyperlinks', 'off'));
    delete(fullfile(process_ID_directory, '*'))
end

end
