%%% Micros Pipeline, Step 2 Check
%%% Following the rereferencing of a microelectrode recording's data
%%% leader node used this function to update progress table with
%%% completion (or error) of this step for each processed recording.
%%% Additionally will add to table signal-to-noise ratios of each recording
%%% before and after referencing, and the channel number selected for reference.

function n02_rereference_signal_check(varargin)
if isempty(varargin)                                                %%% To run manually edit values below
    is_leader_node = true;                                          %%% Always true in this case.
    root_directory = '/path/to/micros_pipeline/parent_directory';   %%% Parent directory.
else                                                                %%% Otherwise this is the order they should be entered into function, following above format
    is_leader_node = varargin{1};
    root_directory = varargin{2};
end

%%% Declare directories
data_directory = fullfile(root_directory, 'micros_database');
check_directory = fullfile(root_directory, 'micros_pipeline/process_files/n02_reref/check');
process_directory = fullfile(root_directory, 'micros_pipeline/process_files/n02_reref/lock_files');
error_directory = fullfile(root_directory, 'micros_pipeline/error_logs');

%%% All nodes wait for last node to finish previous step and create 'ready' file.
ready_yet = false;
while ~ready_yet
    if isfile(fullfile(process_directory, 'ready'))
        ready_yet = true;
    else
        pause(5 + rand)
    end
end

%%% Leader node updates table. Directs stop of the rest of pipeline if there is an error.
start_file = fullfile(check_directory, 'check_start');
done_file = fullfile(check_directory, 'ready');
stop_file = fullfile(check_directory, 'stop');

%%% Leader node deletes file that allows next pipeline step to be started.
%%% Other nodes wait for progress table update to start.
if is_leader_node
    if isfile(done_file)
        delete(done_file)
    end
    if isfile(stop_file)
        delete(stop_file)
    end
else
    check_started = false;
    while ~check_started
        if isfile(start_file)
            check_started = true;
        else
            pause(5 + rand)
        end
    end
end

%%% Leader node updates progress table based on lock, done, and error files generated in previous step.
%%% Other nodes wait for progress table update to end.
if is_leader_node

    file_id = fopen(start_file, 'w'); fclose(file_id);
    fprintf('%s - Node %d performing step 2 check.\n', datestr(now), feature('getpid'));
    
    try
        progress_table_file = fullfile(data_directory, 'progress_table.mat');
        load(progress_table_file, 'progress_table')
        
        lock_files = {dir(fullfile(process_directory, '*_lock')).name}; %list of existing lock_files
        done_files = {dir(fullfile(process_directory, '*_done')).name}; %list of existing lock_files
        error_files = {dir(fullfile(process_directory, 'error_*')).name}; %list of existing lock_files
        
        for ldx = 1:length(lock_files)

            this_lock = lock_files{ldx};
            lock_ID = strrep(this_lock, '_lock', '');
            delimiter = strfind(lock_ID, '_');
            
            subject = lock_ID(1:delimiter(1)-1);
            date = lock_ID(delimiter(1) + 1:delimiter(2) - 1);
            task = lock_ID(delimiter(2) + 1:delimiter(3) - 1);
            folder = lock_ID(delimiter(1) + 1:delimiter(end) - 1);
            bank = lock_ID(end);
            if contains(folder, 'part')       
                part = str2double(folder(end));
            else
                part = 0;
            end
            
            has_subject = contains(progress_table.subject, subject);
            has_date = contains(progress_table.date, date);
            has_task = contains(progress_table.task, task);
            has_part = progress_table.part == part;
            has_bank = contains(progress_table.bank, bank);
            
            row_index = has_subject & has_date & has_task & has_part & has_bank;
            
            is_done = any(contains(done_files, lock_ID));
            had_error = any(contains(error_files, lock_ID));
            
            updated_row = progress_table(row_index, :);
            progress_table(row_index, :) = [];

            if is_done
                
                %%% Get signal-to-noise ratio values before referencing and after referencing
                reref_directory = fullfile(data_directory, subject, folder, 'reref', sprintf('Bank%s', bank));
                load(fullfile(reref_directory, 'SNR_split.mat'), 'SNR_split')
                load(fullfile(reref_directory, 'SNR_reref.mat'), 'SNR_reref')
                
                %%% Get channel number selected as reference
                reref_channel_file = dir(fullfile(data_directory, subject, folder, 'reref', sprintf('Bank%s', bank), 'reref_channel_*'));
                reref_channel_file = {reref_ch_file.name};
                reref_channel_file = reref_ch_file{1};
                reref_channel = regexp(reref_ch_file, '[0-9][0-9][0-9]', 'match');
                reref_channel = str2double(reref_ch{1});
                
                %%% Update values
                updated_row.rerefed = true;
                updated_row.reref_channel = reref_channel;
                updated_row.SNR_split = {SNR_split};
                updated_row.SNR_reref = {SNR_reref};
                
                %%% Delete lock and done files for this entry
                process_done = fullfile(process_directory, done_files{contains(done_files, lock_ID)});
                delete(process_done); 
                delete(fullfile(process_directory, this_lock));
            
            elseif had_error
                
                this_error_directory = fullfile(error_directory, subject, folder, 'reref', sprintf('Bank%s', bank));
                if ~isfolder(this_error_directory)
                    mkdir(this_error_directory);
                end
                
                current_error = fullfile(process_directory, error_files{contains(error_files, lock_ID)});
                error_file = fullfile(this_error_directory, 'error_message.txt');
                movefile(current_error, error_file);
                
                updated_row.has_error = true;

                delete(fullfile(process_directory, this_lock));
            else
                delete(fullfile(process_directory, this_lock));
            end

            progress_table = [progress_table; updated_row];
        end

        progress_table = sortrows(progress_table, {'subject', 'date', 'task', 'part', 'bank'}, 'ascend');
        save(fullfile(data_directory, 'progress_table.mat'), 'progress_table')
        
        %%% In case done before a node finished waiting, pause 10 sec before deleting start file and creating done file.
        pause(10)
        fprintf('%s - Node %d completed step 2 check.\n', datestr(now), feature('getpid'));
        file_id = fopen(done_file, 'w'); fclose(file_id);
        delete(start_file)
        
    catch iteration_error
        %%% In case the error happened before a node finished waiting, pause 10 sec before deleting start file and logging error.
        pause(10)
        fprintf('%s - Node %d errored in step 2 check and quit.\n', datestr(now), feature('getpid'));
        file_id = fopen(stop_file, 'w');
        fprintf(file_id, '%s\n', strcat(datestr(now), ' - Step 2 Check Error:'));
        fprintf(file_id, '%s\n', getReport(iteration_error, 'extended', 'hyperlinks', 'off'));
        fclose(file_id);
        delete(start_file);
        quit
    end
else
    %%% All other nodes will wait for the above to finish.
    %%% Wait periods of 5 sec before rechecking whether the leader node is done updating table.
    check_finished = false;
    while ~check_finished
        pause(5 + rand);
        check_finished = isfile(done_file) | isfile(stop_file);
    end
    if isfile(stop_file)
        fprintf('%s - Node %d quit while waiting for step 2 check completion.\n', datestr(now), feature('getpid'));
        quit
    end
    fprintf('%s - Node %d finished waiting for step 2 check completion.\n', datestr(now), feature('getpid'));
end

end