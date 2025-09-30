%%% Micros Pipeline, Step 3: Filtering Noise
%%% This function will filter line noise from rereferenced signal. 
%%% There are many types of filtering methods, but currently in this
%%% pipeline only utilizing 'zapline' from NoiseTools.
%%% Function will generate a plot with sample signal and power spectral
%%% density for each channel before and after noise filtering.
%%% Using n03_zapline. May add new functions for this step.


function n03_filter_noise(varargin)
if isempty(varargin)                                               %%% To run manually, edit values below.
    is_leader_node = true;                                         %%% Alway true in this case.
    root_directory = '/path/to/micros_pipeline/parent_directory';  %%% Parent directory.
    filtering_method = 'zapline';                                  %%% Currently only 'zapline' from NoiseTools.
else                                                               %%% Otherwise, this is the order they should be entered into function, following above format.
    is_leader_node = varargin{1};
    root_directory = varargin{2};
    filtering_method = varargin{3};
end

%%% Declare directories
data_directory = fullfile(root_directory, 'micros_database');
process_directory = fullfile(root_directory, 'micros_pipeline/process_files/n03_clean');
previous_directory = fullfile(root_directory, 'micros_pipeline/process_files/n02_reref/check');
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

%%% Determine which recordings have been rereferenced, but not cleaned
unfiltered_recordings = progress_table(progress_table.rerefed & ~progress_table.clean, :);

%%% Regular loop to filter unfiltered recordings
%%% Lock files created to avoid multiple node instances from processing the same files.
%%% Not using parallel pool on this step because of high memory usage.

pause(rand * 5) %%% To separate nodes. Each node has unique random number generator seed.

for idx = 1:height(unfiltered_recordings)
    
    this_recording = unfiltered_recordings(idx, :);
    
    subject = this_recording.subject{:};
    date = this_recording.date{:};
    task = this_recording.task{:};
    part = this_recording.part;
    bank = this_recording.bank{:};
    folder = strcat(date, '_', task);
    if part > 0
        folder = strcat(folder, sprintf('_part%d', part));
    end
    
    this_unfiltered_folder = strcat(subject, '_', folder, '_', bank);
    
    lock_file = fullfile(lock_directory, strcat(this_unfiltered_folder, '_lock'));
    done_file = fullfile(lock_directory, strcat(this_unfiltered_folder, '_done'));
    error_file = fullfile(lock_directory, strcat('error_', this_unfiltered_folder));
    
    pause(rand)  %%% In case lock file is searched for simultaneously
    
    if ~isfile(lock_file) 
        file_id = fopen(lock_file, 'w'); fclose(file_id);
        try
            
            switch filtering_method
                case 'zapline'
                    n03_zapline(root_directory, subject, folder, bank);
            end
            
            file_id = fopen(done_file, 'w'); fclose(file_id);
        
        catch iteration_error
            file_id = fopen(error_file, 'w');
            fprintf(file_id, '%s\n', getReport(iteration_error, 'extended', 'hyperlinks', 'off'));
            fclose(file_id);
        end
    end
end

%%% Only last node to finish should set the ready flag for next step to begin.
%%% Number of lock files should match number of done and error files.
%%% Last node also deletes 'reset' file in current process directory.

pause(0.5)

lock_files = dir(fullfile(lock_directory, '*_lock'));
done_files = dir(fullfile(lock_directory, '*_done'));
error_files = dir(fullfile(lock_directory, 'error_*'));

if length(lock_files) == length(done_files) + length(error_files)
    pause(15)
    fprintf('%s - Last node completed step 3.\n', datestr(now));
    ready_flag = fullfile(lock_directory, 'ready');
    file_id = fopen(ready_flag, 'w'); fclose(file_id);
    delete(fullfile(lock_directory, 'reset'));
end

end