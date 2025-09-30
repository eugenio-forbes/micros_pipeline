%%% Micros Pipeline, Step 6: Clustering Signal Spikes
%%% There are many ways in which data can be clustered.
%%% Currently in this pipeline only using Combinato.
%%% More information about Combinato in https://github.com/jniediek/combinato
%%% There are also many different parameters input into the function
%%% like criteria for artifact rejection, noise filtering, recursion, etc.
%%% that ultimately affect the results of spike clustering.
%%% Thus, this function allows for the use of standard parameters (saved in .txt file).
%%% Alternatively, one may loop through different combinations of parameters.
%%% Clustering results could be compared to narrow down set of parameters
%%% that yield more spikes, or avoid over/under-clustering, for example.

function n06_cluster_signal_spikes(varargin)
if isempty(varargin)                                               %%% To run manually, edit values below.
    is_leader_node = true;                                         %%% Alway true in this case.
    root_directory = '/path/to/micros_pipeline/parent_directory';  %%% Parent directory.
    clustering_method = 'combinato_func';                          %%% Only 'combinato_func' available. Open to adding other clustering methods.
    clustering_mode = 'standard';                                  %%% Clustering modes:
                                                                   %%% - 'standard' uses default parameters one time
                                                                   %%% - 'grid-search' loops through 16 parameter sets plus default
else                                                               %%% Otherwise, this is the order they should be entered into function, following above format.
    is_leader_node = varargin{1};
    root_directory = varargin{2};
    clustering_method = varargin{3};
    clustering_mode = varargin{4};
end

warning off

%%% Declare directories
data_directory = fullfile(root_directory, 'micros_database');
process_directory = fullfile(root_directory, 'micros_pipeline/process_files/n06_cluster');
previous_directory = fullfile(root_directory, 'micros_pipeline/process_files/n05_hdf5/check');
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

%%% Determine which recordings have been converted to .hdf5, but not clustered
unclustered_recordings = progress_table(progress_table.hdf5 & ~progress_table.clustered, :);

%%% Parallel pool loop to cluster spikes from all unclustered recordings using n06_combinato_func.
%%% Lock files created to avoid multiple node instances from processing the same files.

pause(rand * 5) %%% To separate nodes. Each node has unique random number generator seed.

parfor idx = 1:height(unclustered_recordings)
    
    this_recording = unclustered_recordings(idx, :);
    
    subject = this_recording.subject{:};
    date = this_recording.date{:};
    task = this_recording.task{:};
    part = this_recording.part;
    bank = this_recording.bank{:};
    folder = strcat(date, '_', task);
    if part > 0
        folder = strcat(folder, sprintf('_part%d', part));
    end
    
    this_unclustered_folder = strcat(subject, '_', folder, '_', bank);
    
    lock_file = fullfile(lock_directory, strcat(this_unclustered_folder, '_lock'));
    done_file = fullfile(lock_directory, strcat(this_unclustered_folder, '_done'));
    error_file = fullfile(lock_directory, strcat('error_', this_unclustered_folder));
    
    pause(rand)  %%% In case lock file is searched for simultaneously
    
    if ~isfile(lock_file) 
        file_id = fopen(lock_file, 'w'); fclose(file_id);
        try
        
            switch clustering_method
                case 'combinato_func'
                    [error_flag, error_message] = n06_combinato_func(root_directory, subject, folder, bank, clustering_mode);
            end
            
            %%% System commands are used in n06_combinato_func to run python-based Combinato
            %%% Catching errors based on system output
            
            if error_flag
                file_id = fopen(error_file, 'w');
                fprintf(file_id, '%s\n', error_message);
                fclose(file_id);
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

%%% Delete the parallel pool until step 8 is completed.

%%% Only last node to finish should set the ready flag for next step to begin.
%%% Number of lock files should match number of done and error files.
%%% Last node also deletes 'reset' file in current process directory.

pause(0.5)

lock_files = dir(fullfile(lock_directory, '*_lock'));
done_files = dir(fullfile(lock_directory, '*_done'));
error_files = dir(fullfile(lock_directory, 'error_*'));

if length(lock_files) == length(done_files) + length(error_files)
    pause(15)
    fprintf('%s - Last node completed step 6.\n', datestr(now));
    ready_flag = fullfile(lock_directory, 'ready');
    file_id = fopen(ready_flag, 'w'); fclose(file_id);
    delete(fullfile(lock_directory, 'reset'));
end

end
