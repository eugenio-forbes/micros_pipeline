%%% Micros Pipeline, Step 11 Check
%%% Following the alignment of behavioral data to a microelectrode recording's data,
%%% leader node uses this function to update progress table with
%%% completion (or error) of this step for each processed recording.
%%% Additionally will add to table R-squared and max deviation of results of alignment.

function n11_align_behavioral_data_check(varargin)
if isempty(varargin)                                                %%% To run manually edit values below
    is_leader_node = true;                                          %%% Always true in this case.
    root_directory = '/path/to/micros_pipeline/parent_directory';   %%% Parent directory.
else                                                                %%% Otherwise this is the order they should be entered into function, following above format
    is_leader_node = varargin{1};
    root_directory = varargin{2};
end

%%% Declare directories
data_directory = fullfile(root_directory, 'micros_database');
check_directory = fullfile(root_directory, 'micros_pipeline/process_files/n11_align/check');
process_directory = fullfile(root_directory, 'micros_pipeline/process_files/n11_align/lock_files');
completion_directory = fullfile(root_directory, 'micros_pipeline/completion_reports');
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
    fprintf('%s - Node %d performing step 11 check.\n', datestr(now), feature('getpid'));
    
    try
        progress_table_file = fullfile(data_directory, 'progress_table.mat');
        load(progress_table_file, 'progress_table')
        
        lock_files = {dir(fullfile(process_directory, '*_lock')).name};
        done_files = {dir(fullfile(process_directory, '*_done')).name};
        error_files = {dir(fullfile(process_directory, 'error_*')).name};
        
        for ldx = 1:length(lock_files)

            this_lock = lock_files{ldx};
            lock_ID = strrep(this_lock, '_lock', '');
            delimiter = strfind(lock_ID, '_');
            
            subject = lock_ID(1:delimiter(1)-1);
            date = lock_ID(delimiter(1)+1:delimiter(2)-1);
            folder = lock_ID(delimiter(1)+1:end);
            if contains(folder, 'part')
                task = lock_ID(delimiter(2)+1:delimiter(3)-1);
                part = str2double(folder(end));
            else
                task = lock_ID(delimiter(2)+1:end);
                part = 0;
            end
            
            has_subject = contains(progress_table.subject, subject);
            has_date = contains(progress_table.date, date);
            has_task = contains(progress_table.task, task);
            has_part = progress_table.part == part;
            
            row_indices = has_subject & has_date & has_task & has_part;
            
            is_done = any(contains(done_files, lock_ID));
            had_error = any(contains(error_files, lock_ID));
            
            updated_rows = progress_table(row_indices, :);
            progress_table(row_indices, :) = [];
            n_microelectrodes = height(updated_rows);

            if is_done
            
                %%% Get alignment report stats
                align_directory = fullfile(data_directory, subject, folder, 'alignment');
                report_file = fullfile(align_directory, 'report_stats.mat');
                load(report_file, 'report_stats')
                
                %%% Update values
                updated_rows.r_squared = repmat(report_stats.r_squared, n_microelectrodes, 1);
                updated_rows.max_deviation = repmat(report_stats.max_deviation, n_microelectrodes, 1);
                updated_rows.aligned = true(n_microelectrodes, 1);
                
                process_done = fullfile(process_directory, done_files{contains(done_files, lock_ID)});
                delete(process_done);
                delete(fullfile(process_directory, this_lock));
                
            elseif had_error
            
                this_error_directory = fullfile(error_directory, subject, folder, 'align');
                if ~isfolder(this_error_directory)
                    mkdir(this_error_directory);
                end
                
                current_error = fullfile(process_directory, error_files{contains(error_files, lock_ID)});
                error_file = fullfile(this_error_directory, 'error_message.txt');
                movefile(current_error, error_file);
                
                updated_rows.has_error = true(n_microelectrodes, 1);

                delete(fullfile(process_directory, this_lock));
            else
                delete(fullfile(process_directory, this_lock));
            end

            progress_table = [progress_table; updated_rows];
        end

        progress_table = sortrows(progress_table, {'subject', 'date', 'task', 'part', 'bank'}, 'ascend');
        save(fullfile(data_directory, 'progress_table.mat'), 'progress_table')
        
        %%% In case done before a node finished waiting, pause 10 sec before deleting start file and creating done file.
        pause(10)
        fprintf('%s - Node %d completed step 10 check.\n', datestr(now), feature('getpid'));
        file_id = fopen(done_file, 'w'); fclose(file_id);
        delete(start_file)
        
    catch iteration_error
        %%% In case the error happened before a node finished waiting, pause 10 sec before deleting start file.
        pause(10)
        fprintf('%s - Node %d errored in step 10 check and quit.\n', datestr(now), feature('getpid'));
        file_id = fopen(stop_file, 'w');
        fprintf(file_id, '%s\n', strcat(datestr(now), ' - Step 10 Check Error:'));
        fprintf(file_id, '%s\n', getReport(iteration_error, 'extended', 'hyperlinks', 'off'));
        fclose(file_id);
        delete(start_file);
        quit
    end
else
    %%% All other nodes will wait for the above to finish.
    %%% Wait periods of 5 sec before rechecking whether the chosen node is done updating table.
    check_finished = false;
    while ~check_finished
        pause(5 + rand);
        check_finished = isfile(done_file) | isfile(stop_file);
    end
    if isfile(stop_file)
        fprintf('%s - Node %d quit while waiting for step 11 check completion.\n', datestr(now), feature('getpid'));
        quit
    end
    fprintf('%s - Node %d finished waiting for step 11 check completion.\n', datestr(now), feature('getpid'));
end

end