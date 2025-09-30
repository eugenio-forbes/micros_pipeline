%%% Micros Pipeline, Step 11: Aligning Behavioral Data
%%% For behavioral tasks that delivered sync pulses to EEG recording,
%%% behavioral data computer time is aligned to EEG recording time,
%%% so that spike data analyses related to behavior can be performed.

function n11_align_behavioral_data(varargin)
if isempty(varargin)                                               %%% To run manually, edit values below
    is_leader_node = true;                                         %%% Alway true in this case
    n_workers_align = 10;                                          %%% Number of workers aligning behavioral data in parallel pool.
    root_directory = '/path/to/micros_pipeline/parent_directory';  %%% Parent directory.
else                                                               %%% Otherwise, this is the order they should be entered into function, following above format
    is_leader_node = varargin{1};
    n_workers_align = varargin{2};
    root_directory = varargin{3};
end

%%% Declare directories
data_directory = fullfile(root_directory, 'micros_database');
process_directory = fullfile(root_directory, 'micros_pipeline/process_files/n11_align');
previous_directory = fullfile(root_directory, 'micros_pipeline/process_files/n10_modal/check');
lock_directory = fullfile(process_directory, 'locks');

%%% List tasks with behavioral events files
event_tasks = {'AR', 'AR-scopolamine', 'AR-stim', 'FR', 'FR-scopolamine', 'SR'};

%%% All nodes wait for leader node to complete updating progress table after checking processing of previous step.
ready_yet = false;
while ~ready_yet
    if isfile(fullfile(previous_directory, 'ready'))
        ready_yet = true;
    elseif isfile(fullfile(previous_directory, 'stop'))
        quit
    else
        pause(5 + rand)
    end
end

%%% Leader node resets lock files in process directory, if any remaining.
%%% Other nodes wait until the previously existing lock files have been reset.
reset_file = fullfile(lock_directory, 'reset');
if is_leader_node
    lock_files = dir(lock_directory);
    lock_files = lock_files(~contains({lock_files.name}, {'.', '..'}));
    if ~isempty(lock_files)
        delete(fullfile(lock_directory, '*'));
    end
    file_id = fopen(reset_file, 'w'); fclose(file_id);
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

%%% Determine which recordings with events have not been aligned. Can only align those with sync files. 
unaligned_recordings = progress_table(progress_table.events & progress_table.has_sync_pulses & ~progress_table.aligned, :);
[~, unique_sessions, ~] = unique(fullfile(unaligned_recordings.subject, unaligned_recordings.date, unaligned_recordings.task, num2str(unaligned_recordings.part)));
unaligned_recordings = unaligned_recordings(unique_sessions, :);

%%% Parallel pool loop to rereference all unrereferenced recordings using n02_rereference_func.
%%% Lock files created to avoid multiple node instances from processing the same files.
parpool('local', n_workers_align, 'IdleTimeout', 1440)

pause(rand * 5) %%% To separate nodes. Each node has unique random number generator seed.

parfor idx = 1:height(unaligned_recordings)
    
    this_recording = unaligned_recordings(idx, :);
    
    subject = this_recording.subject{:};
    date = this_recording.date{:};
    task = this_recording.task{:};
    part = this_recording.part;
    has_events = this_recording.has_events;
    sync_channel = this_recording.sync_channel;
    folder = strcat(date, '_', task);
    if part > 0
        folder = strcat(folder, sprintf('_part%d', part));
    end
    
    unaligned_folder = strcat(subject, '_', folder);
    
    lock_file = fullfile(lock_directory, strcat(unaligned_folder, '_lock'));
    done_file = fullfile(lock_directory, strcat(unaligned_folder, '_done'));
    error_file = fullfile(lock_directory, strcat('error_', unaligned_folder));
    
    pause(rand)  %%% In case lock file is searched for simultaneously
    
    if ~isfile(lock_file) 
        file_id = fopen(lock_file, 'w'); fclose(file_id);
        try
        
            if ismember(task, event_tasks)
            
                if has_events
                
                    [error_flag, error_message] = n10_align_func(root_directory, subject, folder, sync_channel);
                    
                    if error_flag
                        file_id = fopen(error_file, 'w'); 
                        fprintf(file_id, error_message);
                        fclose(file_id);
                    else
                        file_id = fopen(done_file, 'w'); fclose(file_id);
                    end
                    
                else
                    file_id = fopen(error_file, 'w'); 
                    fprintf(file_id, 'Missing events file.');
                    fclose(file_id);
                end
                
            else
                file_id = fopen(done_file, 'w'); fclose(file_id);
            end
            
        catch iteration_error
            file_id = fopen(error_file, 'w');
            fprintf(file_id, '%s\n', getReport(iteration_error, 'extended', 'hyperlinks', 'off'));
            fclose(file_id);
        end
    end
end

%%% Delete parallel pool.
delete(gcp('nocreate'))

%%% Only last node to finish should set the ready flag for next step to begin.
%%% Number of lock files should match number of done and error files.
%%% Last node also deletes 'reset' file in current process directory.

pause(0.5)

lock_files = dir(fullfile(lock_directory, '*_lock'));
done_files = dir(fullfile(lock_directory, '*_done'));
error_files = dir(fullfile(lock_directory, 'error_*'));

if length(lock_files) == length(done_files) + length(error_files)
    pause(15)
    fprintf('%s - Last node completed step 11.\n', datestr(now));
    ready_flag = fullfile(lock_directory, 'ready');
    file_id = fopen(ready_flag, 'w'); fclose(file_id);
    delete(fullfile(lock_directory, 'reset'));
end

end