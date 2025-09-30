%%% Micros Pipeline, Step 1: Splitting Microelectrode Recordings into Separate Files
%%% This function will gather a list of all raw recordings (.ns6 files).
%%% Based on a saved progress table it will determine which recordings have not been split
%%% Each microelectrode carries data for 8 channels. A recording may contain data from more
%%% than one microelectrode, and additionally it may contain data from a DC channel used
%%% for delivering pulses that would align the recording to the behavioral data.
%%% Thus, the data of each recording will be split into individual files for each channel,
%%% using n01_split_ns6.m

function n01_split_recording_file(varargin)
if isempty(varargin)                                               %%% To run manually edit values below.
    time_threshold = 0;                                            %%% In this case alway 0.
    is_leader_node = true;                                         %%% Always true in this case.
    n_workers_split = 2;                                           %%% Number of workers in parallel pool splitting recordings.
    root_directory = '/path/to/micros_pipeline/parent_directory';  %%% Parent directory.
else                                                               %%% Otherwise this is the order they should be entered into function, following above format.
    time_threshold = varargin{1};
    is_leader_node = varargin{2};
    n_workers_split = varargin{3};
    root_directory = varargin{4};
end

%%% Declare directories
data_directory = fullfile(root_directory, 'micros_database');
process_directory = fullfile(root_directory, 'micros_pipeline/process_files/n01_split');
lock_directory = fullfile(process_directory, 'locks');

%%% Leader node resets lock files in process directory, if any remaining.
%%% Other nodes wait until the previously existing lock files have been reset.
reset_file = fullfile(lock_directory, 'reset');
if is_leader_node
    lock_files = dir(lock_directory);
    lock_files = lock_files(~contains({lock_files.name}, {'.', '..'}));
    if ~isempty(lock_files)
        delete(fullfile(lock_directory, '*'));
    end
    file_id = fopen(reset_file, 'w');
    fclose(file_id);
else
    reset_completed = false;
    while ~reset_completed
        if isfile(reset_file)
            reset_completed = true;
        else
            pause(5 + rand)
        end
    end
end

%%% Load table with information about progress in processing micro recordings
progress_table_file = fullfile(data_directory, 'progress_table.mat');
load(progress_table_file, 'progress_table')

%%% Gather information of all .ns6 files in the database (including new ones added manually)
ns6_paths = dir(fullfile(data_directory, '*/*/raw/*.ns6'));
ns6_paths = unique(fullfile({ns6_paths.folder}))';
ns6_folders = cell(length(ns6_paths), 1);
ns6_subjects = cell(length(ns6_paths), 1);
for idx = 1:length(ns6_paths)
    this_ns6_path = ns6_paths{idx};
    delimiter = strfind(this_ns6_path,'/');
    folder_name = this_ns6_path(delimiter(end - 1) + 1:delimiter(end) - 1);
    subject = this_ns6_path(delimiter(end - 2) + 1:delimiter(end - 1) - 1);
    ns6_subjects{idx} = subject;
    ns6_folders{idx} = folder_name;
end

%%% Determine from the progress table which files have already been split
split_recordings = progress_table(progress_table.split, :);
split_recording_paths = cell(height(split_recordings), 1);
for idx = 1:height(split_recordings)
    subject = split_recordings(idx, :).subject;
    task = split_recordings(idx, :).task;
    date = split_recordings(idx, :).date;
    part = split_recordings(idx, :).part;
    folder_name = strcat(date, '_', task);
    
    if part > 0   %%% For test days with multiple parts for a session or multiple sessions
        folder_name = strcat(folder, sprintf('_part%d', part));
    end
    
    split_recording_paths{idx} = fullfile(data_directory, subject{:}, folder_name{:}, 'raw');
end

%%% Determine which recordings have not been split
unsplit_recordings = ~ismember(ns6_paths, split_recording_paths);
unsplit_subjects = ns6_subjects(unsplit_recordings);
unsplit_folders = ns6_folders(unsplit_recordings);

%%% Parallel pool loop to split all unsplit recordings using n01_split_ns6.
%%% Lock files created to avoid multiple node instances from processing the same files.

parpool('local', n_workers_split, 'IdleTimeout', 1440)

pause(rand * 5) %%% To separate nodes. Each node has unique random number generator seed.

parfor idx = 1:length(unsplit_folders)
    
    this_subject = unsplit_subjects{idx};
    this_folder = unsplit_folders{idx};
    this_unsplit_folder = strcat(this_subject, '_', this_folder);
    
    lock_file = fullfile(lock_directory, strcat(this_unsplit_folder, '_lock'));
    done_file = fullfile(lock_directory, strcat(this_unsplit_folder, '_done'));
    error_file = fullfile(lock_directory, strcat('error_', this_unsplit_folder));
    
    pause(rand) %%% In case lock file is searched for simultaneously
    
    if ~isfile(lock_file)
        file_id = fopen(lock_file, 'w'); fclose(file_id);
        try
        
            n01_split_ns6(root_directory, this_subject, this_folder);
            
            file_id = fopen(done_file, 'w'); fclose(file_id);
            
        catch iteration_error
            file_id = fopen(error_file, 'w');
            fprintf(file_id, '%s\n', getReport(iteration_error, 'extended', 'hyperlinks', 'off'));
            fclose(file_id);
        end
    end
end

%%% Deleting parallel pool after each step to start new ones with different number of workers.
delete(gcp('nocreate'))

%%% For processing steps, errors are not detrimental to the functioning of the pipeline. 
%%% Errors will be logged into table during check steps.
%%% Full pipeline stop should only occur for errors during check steps, where the table is updated.

%%% Only last node to finish should set the ready flag for next step to
%%% begin. So we check for number of lock files to match number of done and
%%% error files. Then we create a ready flag and delete reset lock.

pause(0.5)
lock_files = dir(fullfile(lock_directory, '*_lock'));
done_files = dir(fullfile(lock_directory, '*_done'));
error_files = dir(fullfile(lock_directory, 'error_*'));

if length(lock_files) == length(done_files) + length(error_files)
    %%% Since this is the first step, pause for a time equivalent to the
    %%% time_threshold set in initialization to allow for nodes starting
    %%% within this time interval to catch up before the reset is deleted,
    %%% so they can be synced for next step.
    pause(time_threshold + 15)
    fprintf('%s - Last node completed step 1.\n', datestr(now));
    ready_flag = fullfile(lock_directory, 'ready');
    file_id = fopen(ready_flag, 'w'); 
    fclose(file_id);
    delete(fullfile(lock_directory, 'reset'));
end

end
