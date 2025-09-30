%%% Micros Pipeline, Step 4: Rescaling Signal
%%% Clustering method 'combinato' assumes that input signal is in uV.
%%% EEG data from Blackrock Neurotech saved in int16 format needs to be
%%% divided by a factor of 4 in order to be in uV. This leads to a loss
%%% of resolution when keeping data in int16 format. Saving a copy
%%% to be used for spike clustering by combinato only. Can use this
%%% copy using combinato's Matlab data input option or .hdf5 file.
%%% Combinato works best when using .hdf5 files.

function n04_rescale_signal(varargin)
if isempty(varargin)                                               %%% To run manually, edit values below.
    is_leader_node = true;                                         %%% Alway true in this case.
    n_workers_rescale = 2;                                         %%% Number of workers performing rescaling in parallel pool.
    root_directory = '/path/to/micros_pipeline/parent_directory';  %%% Parent directory.
else                                                               %%% Otherwise, this is the order they should be entered into function, following above format.
    is_leader_node = varargin{1};
    n_workers_rescale = varargin{2};
    root_directory = varargin{3};
end

%%% Declare directories
data_directory = fullfile(root_directory, 'micros_database');
process_directory = fullfile(root_directory, 'micros_pipeline/process_files/n04_rescale');
previous_directory = fullfile(root_directory, 'micros_pipeline/process_files/n03_clean/check');
lock_directory = fullfile(process_directory, 'locks');

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

%%% Determine which recordings have been cleaned, but not rescaled
unrescaled_recordings = progress_table(progress_table.clean & ~progress_table.rescaled, :);

%%% Parallel pool loop to rescale all unrescaled recordings using n04_rescale_func.
%%% Lock files created to avoid multiple node instances from processing the same files.
parpool('local', n_workers_rescale, 'IdleTimeout', 1440)

pause(rand * 5) %%% To separate nodes. Each node has unique random number generator seed.

parfor idx = 1:height(unrescaled_recordings)
    
    this_recording = unrescaled_recordings(idx, :);
    
    subject = this_recording.subject{:};
    date = this_recording.date{:};
    task = this_recording.task{:};
    part = this_recording.part;
    bank = this_recording.bank{:};
    folder = strcat(date, '_', task);
    if part > 0
        folder = strcat(folder, sprintf('_part%d', part));
    end
    
    this_unrescaled_folder = strcat(subject, '_', folder, '_', bank);
    
    lock_file = fullfile(lock_directory, strcat(this_unrescaled_folder, '_lock'));
    done_file = fullfile(lock_directory, strcat(this_unrescaled_folder, '_done'));
    error_file = fullfile(lock_directory, strcat('error_', this_unrescaled_folder));
    
    pause(rand)  %%% In case lock file is searched for simultaneously
    
    if ~isfile(lock_file) 
        file_id = fopen(lock_file, 'w'); fclose(file_id);
        try
            
            n04_rescale_func(root_directory, subject, folder, bank);
            
            file_id = fopen(done_file, 'w'); fclose(file_id);
            
        catch iteration_error
            file_id = fopen(error_file, 'w');
            fprintf(file_id, '%s\n', getReport(iteration_error, 'extended', 'hyperlinks', 'off'));
            fclose(file_id);
        end
    end
end

%%% Delete parallel pool so that a new one with different number of workers can be started.
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
    fprintf('%s - Last node completed step 4.\n', datestr(now));
    ready_flag = fullfile(lock_directory, 'ready');
    file_id = fopen(ready_flag, 'w'); fclose(file_id);
    delete(fullfile(lock_directory, 'reset'));
end

end