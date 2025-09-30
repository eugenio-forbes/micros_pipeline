%%% Micros Pipeline, Step 1 Check
%%% Following the splitting of all raw recordings into the individual files
%%% for each of the channels, this function will gather information for each
%%% newly split recording to add to the progress table. Channels for each
%%% recording are assigned to a 'bank' corresponding to microelectrode in
%%% Blackrock Neurotech hardware. Function gathers length of recordings,
%%% whether a channel contains sync pulses, whether there are behavioral events
%%% associated to recording, location of microelectrodes, whether subject participated 
%%% in UPenn experiments. 

function n01_split_recording_file_check(varargin)
if isempty(varargin)                                                %%% To run manually edit values below
    is_leader_node = true;                                          %%% Always true in this case.
    root_directory = '/path/to/micros_pipeline/parent_directory';   %%% Parent directory.
else                                                                %%% Otherwise this is the order they should be entered into function, following above format
    is_leader_node = varargin{1};
    root_directory = varargin{2};
end

%%% Declare directories
data_directory = fullfile(root_directory, 'micros_database');
check_directory = fullfile(root_directory, 'micros_pipeline/process_files/n01_split/check');
process_directory = fullfile(root_directory, 'micros_pipeline/process_files/n01_split/locks');
error_directory = fullfile(root_directory, 'micros_pipeline/error_logs');

%%% Assuming a maximum of 16 micro depth electrodes would be implanted.
%%% So far, there has been 1-3 implanted. Could achieve maximum of 32 if we had the equipment.
bank_names = {'A'; 'B'; 'C'; 'D'; 'E'; 'F'; 'G'; 'H'; 'I'; 'J'; 'K'; 'L'; 'M'; 'N'; 'O'; 'P'};

%%% List tasks that have behavioral events to align
event_tasks = {'AR', 'AR-scopolamine', 'AR-stim', 'FR', 'FR-scopolamine', 'SR'};

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
    fprintf('%s - Node %d performing step 1 check.\n', datestr(now), feature('getpid'));
    
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
            
            subject = lock_ID(1:delimiter(1) - 1);
            date = lock_ID(delimiter(1) + 1:delimiter(2) - 1);
            folder = lock_ID(delimiter(1) + 1:end);
            if contains(folder, 'part')
                task = lock_ID(delimiter(2) + 1:delimiter(3) - 1);            
                part = str2double(lock_ID(end));
            else
                task = lock_ID(delimiter(2) + 1:end);
                part = 0;
            end
            
            has_subject = contains(progress_table.subject, subject);
            has_date = contains(progress_table.date, date);
            has_task = contains(progress_table.task, task);
            has_part = progress_table.part == part;
            
            row_indices = has_subject & has_date & has_task & has_part;
            
            is_done = any(contains(done_files, lock_ID));
            had_error = any(contains(error_files, lock_ID));
            
            %%% Create new rows or update new ones based on this information
            if sum(row_indices) == 0
                new_entry = progress_table(1, :); %%% First row is a template. Reset to add new entry.
                new_entry.has_error = false;      %%% Template marked with error so it can be excluded.
                new_entry.penn = 0; 
                new_entry.penn_ID = {'O'};        %%% Change manually if not found by code
                new_entry.subject = {subject};
                new_entry.date = {date};
                new_entry.task = {task};
                new_entry.part = part;
                new_entry.notes = {struct('note', {'new_entry'})};
                n_microelectrodes = determine_n_microelectrodes(data_directory, subject, folder);
                
                if had_error || n_microelectrodes == 0
                    new_rows = new_entry;
                    new_rows.errors = {struct('error', {'Not split'})};
                    new_rows.n_microelectrodes = 0;
                else
                    new_rows = repmat(new_entry, n_microelectrodes, 1);
                    new_rows.bank = bank_names(1:n_microelectrodes);
                    new_rows.n_microelectrodes = repmat(n_microelectrodes, n_microelectrodes, 1);
                    new_rows.notes = repelem({struct('note', {''})}, n_microelectrodes, 1);
                end
                
            else
                old_rows = progress_table(row_indices, :);
                progress_table(row_indices, :) = [];
                
                if sum(row_indices) == 1
                    table_notes = old_rows.notes{1};
                    table_notes = {table_notes.note};
                    if any(contains(table_notes, 'new_entry'))         %%% New entry that previously had error
                        n_microelectrodes = determine_n_microelectrodes(data_directory, subject, folder);
                        if had_error || n_microelectrodes == 0         %%% Keep it the same if there was error again
                            new_rows = old_rows;
                        else                                           %%% Delete error from row if it worked this time.
                            new_rows = repmat(old_rows, n_microelectrodes, 1);
                            new_rows.bank = bank_names(1:n_microelectrodes);
                            new_rows.n_microelectrodes = repmat(n_microelectrodes, n_microelectrodes, 1);
                            new_rows.notes = repelem({struct('note', {''})}, n_microelectrodes, 1);
                            new_rows.errors = repelem({struct('error', {''})}, n_microelectrodes, 1);
                        end
                        
                    else                                               %%% Recording with only one bank from previous table entry with no errors.
                        new_rows = old_rows;
                    end
                    
                else %%% Assuming it was a previous entry with more than one bank that had already been entered in table with no errors.
                    new_rows = old_rows;
                end
            end
            
            %%% Now that new entries have been created (or old entries identified) information can be added to new rows and the table can be updated.
            if is_done  %%% Only change values of new rows if split is successful
                
                %%% Get all information about .ns6 recording split
                [file_length, sampling_rate, n_microelectrodes, labels, channels, locations, has_sync_pulses, sync_channel,...
                    sync_file, is_digital, has_events] = get_splitting_info(data_directory, subject, folder);
                [is_penn_subject, penn_ID] = determine_if_penn_subject(subject, n_microelectrodes);
                
                %%% Update values
                new_rows.split = true(n_microelectrodes, 1);                                 %%% Logical indicating process completed
                new_rows.bank = bank_names(1:n_microelectrodes);                             %%% Letter depending on the order of the channels
                new_rows.label = labels;                                                     %%% Depth electrode label
                new_rows.channels = channels;                                                %%% Blackrock channels where micros were recorded
                new_rows.location = locations;                                               %%% As identified by macro channel 1 o micro macro in depth_el_info
                new_rows.events = repmat(ismember(task, event_tasks), n_microelectrodes, 1); %%% Whether there should be events
                new_rows.has_events = has_events;                                            %%% Whether there are events files in folder.
                new_rows.has_sync_pulses = has_sync_pulses;                                  %%% Split file has ainp channel with sync pulses?
                new_rows.sync_channel = sync_channel;                                        %%% If it does what channel number?
                new_rows.sync_file = sync_file;                                              %%% NS6 or NS3?
                new_rows.threshold = repmat(5, n_microelectrodes, 1);                        %%% In case threshold is eventually reduced in future attempts.
                new_rows.digital_headstage = is_digital;                                     %%% Based on channel numbers. Will have to verify for parkland subjects.
                new_rows.penn = is_penn_subject;                                             %%% Based on the existance of behavioral or raw .json files containing R0000T pattern.
                new_rows.penn_ID = penn_ID;                                                  %%% ID found
                new_rows.file_length = file_length;                                          %%% To speed up processing of files
                new_rows.sampling_rate = sampling_rate;                                      %%% Sampling rate
                
                %%% Delete lock and done files for this entry
                process_done = fullfile(process_directory, done_files{contains(done_files, lock_ID)});
                delete(process_done); 
                delete(fullfile(process_directory, this_lock));
            
            elseif had_error %%% Move and rename error file to error logs folder
                this_error_directory = fullfile(error_directory, subject, folder, 'split');
                if ~isfolder(this_error_directory)
                    mkdir(this_error_directory);
                end
                current_error = fullfile(process_directory, error_files{contains(error_files, lock_ID)});
                error_file = fullfile(this_error_directory, 'error_message.txt');
                movefile(current_error, error_file);
                new_rows.has_error = true(height(new_rows), 1);
                %%% Delete lock
                delete(fullfile(process_directory, this_lock));
                
            else %%% The execution of the pipeline did not complete, only creating lock file, so will only delete lock
                delete(fullfile(process_directory, this_lock));
            end
            
            %%% Add new updated rows to table
            progress_table = [progress_table; new_rows];
        end
        
        progress_table = sortrows(progress_table, {'subject', 'date', 'task', 'part', 'bank'}, 'ascend');
        save(fullfile(data_directory, 'progress_table.mat'), 'progress_table')
       
        %%% In case done before a node finished waiting, pause 10 sec before deleting start file and creating done file.
        pause(10)
        fprintf('%s - Node %d completed step 1 check.\n', datestr(now), feature('getpid'));
        file_id = fopen(done_file, 'w'); fclose(file_id);
        delete(start_file)
        
    catch iteration_error
        %%% In case the error happened before a node finished waiting, pause 10 sec before deleting start file.
        pause(10)
        fprintf('%s - Node %d errored in step 1 check and quit.\n', datestr(now), feature('getpid'));
        file_id = fopen(stop_file, 'w');
        fprintf(file_id, '%s\n', strcat(datestr(now), ' - Step 1 Check Error:'));
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
        fprintf('%s - Node %d quit while waiting for step 1 check completion.\n', datestr(now), feature('getpid'));
        quit
    end
    fprintf('%s - Node %d finished waiting for step 1 check completion.\n', datestr(now), feature('getpid'));
end

end


%%% Function to determine number of microelectrodes in the recording (8 channels per microelectrode)
function n_microelectrodes = determine_n_microelectrodes(data_directory, subject, folder)
split_directory = fullfile(data_directory, subject, folder, 'split');
file_list = dir(split_directory);
file_list = file_list(~[file_list.isdir]);
file_list = file_list(~contains({file_list.name}, {'.txt', '.ns6', 'time'}));
n_microelectrodes = idivide(int16(length(file_list)), int16(8));
end


%%% Function to gather information about recordings: file length, sampling rate, number of microelectrodes,
%%% labels, channel numbers, microelectrode brain location, whether there is a DC channel for sync pulses in recording,
%%% whether the hardware used is analog or digital, and whether there are behavioral events files associated to recording.
function [file_length, sampling_rate, n_microelectrodes, labels, channels, locations, has_sync_pulses, sync_channel,...
    sync_file, is_digital, has_events] = get_splitting_info(data_directory, subject, folder)
    
split_directory = fullfile(data_directory, subject, folder, 'split');
file_list = dir(split_directory);
files = fullfile({file_list.folder}, {file_list.name});
jacksheet_file = files{contains({file_list.name}, 'jacksheet')};
length_file = files{contains({file_list.name}, 'file_length')};
files = file_list(~contains({file_list.name}, {'.txt', '.ns6', 'time', 'length'}) & ~ismember({file_list.name}, {'.', '..'}));
files = {files.name};
channel_numbers = regexp(files, '[0-9][0-9][0-9]', 'match');
channel_numbers = cellfun(@(x) str2double(x), channel_numbers);
n_microelectrodes = idivide(int16(length(files)), int16(8));
recording_measures = load(length_file);
file_length = recording_measures.file_length(1);
sampling_rate = recording_measures.sampling_rate(1);
file_length = repmat(file_length, n_microelectrodes, 1);
sampling_rate = repmat(sampling_rate, n_microelectrodes, 1);

%%% Initialize variables for depth electrode data and set default to the case where the data is either not in the right format or not found.
labels = repelem({'Enter manually'}, n_microelectrodes, 1);
locations = repelem({'Enter manually'}, n_microelectrodes, 1);
channels = cell(n_microelectrodes, 1);
has_sync_pulses = false(n_microelectrodes, 1);
sync_channel = zeros(n_microelectrodes, 1);
sync_file = repelem({'NS0'}, n_microelectrodes, 1);
is_digital = false(n_microelectrodes, 1);

%%% Get labels and channels from jacksheet, get location from clinical depth electrode information text file.
if isfile(jacksheet_file)
    jacksheet = readtable(jacksheet_file, delimitedTextImportOptions('Delimiter', ' '), 'ReadVariableNames', false);
    jacksheet_labels = unique(jacksheet{:, 2}, 'stable');
    
    if any(contains(jacksheet_labels, 'ainp')) %%% Analog input channel corresponding to DC channels
        has_sync_pulses = true(n_microelectrodes, 1);
        sync_channel = repmat(str2double(jacksheet{contains(jacksheet{:, 2}, 'ainp'), 1}), n_microelectrodes, 1);
        channel_numbers = channel_numbers(channel_numbers~=sync_channel(1));
        sync_file = repelem({'NS6'}, n_microelectrodes, 1);
    end
    
    jacksheet_labels = jacksheet_labels(~contains(jacksheet_labels, 'ainp'));
    labels = cellfun(@(x) regexprep(x, '\d', ''), jacksheet_labels, 'UniformOutput', false);
    labels = unique(labels, 'stable');
    
    if length(labels) == n_microelectrodes
    
        for bdx = 1:n_microelectrodes
            this_label = labels{bdx};
            microelectrode_channel_numbers = cellfun(@(x) str2double(x), jacksheet{contains(jacksheet{:, 2}, this_label), 1}, 'UniformOutput', true);
            min_channel_number = min(microelectrode_channel_numbers);
            max_channel_number = max(microelectrode_channel_numbers);
            if min_channel_number > 128 %%% All analog recordings below 128. All digital recordings above 128 except Parkland subjects. Verify.
                is_digital(bdx) = true;
            end
            channels{bdx} = [min_channel_number, max_channel_number];
            locations{bdx} = get_microelectrode_location(subject, this_label);
        end
        
    elseif length(channel_numbers) == (n_microelectrodes * 8)
        labels = repelem({'OO'}, n_microelectrodes, 1);
        
        for bdx = 1:n_microelectrodes
            channels{bdx} = [channel_numbers((bdx * 8) - 7), chan_numbers(bdx * 8)];
        end
        
    else
        labels = repelem({'OO'}, n_microelectrodes, 1);
        channels = repelem({[0, 0]}, n_microelectrodes, 1);
    end
    
else
    labels = repelem({'OO'}, n_microelectrodes, 1);
    channels = repelem({[0, 0]}, n_microelectrodes, 1);
end

%%% Check whether the events file is found in the same session folder as
%%% the .ns6. Could potentially identify other locations based on creation
%%% date and which part is being processed. But there could be two parts
%%% with micro recordings and one part without.
events_file = fullfile(data_directory, subject, folder, 'behavioral', 'events.mat');
has_events = repmat(isfile(events_file), n_microelectrodes, 1);

end


%%% Function to determine microelectrode location assuming that it will be the same as location
%%% of macro-micro depth electrode's channel 1 (tip of macro electrode).
function location = get_microelectrode_location(subject, label)
label = strcat(label, '1');
location = 'Enter manually';
subject_directory = sprintf('/path/to/subject_files/%s/documents', subject);
depth_electrode_info_file = dir(fullfile(subject_directory, '*depth*'));
depth_electrode_info_file = depth_electrode_info_file(~contains({depth_electrode_info_file.name}, '._'));

if ~isempty(depth_file)
    depth_electrode_info_file = depth_electrode_info_file(1);
    depth_electrode_info_file = fullfile({depth_electrode_info_file.folder}, {depth_electrode_info_file.name});
    depth_electrode_info_file = depth_electrode_info_file{1};
    depth_electrode_info = readtable(depth_electrode_info_file, delimitedTextImportOptions('Delimiter', ' '), 'ReadVariableNames', false);
    index = find(contains(depth_electrode_info{:, 2}, label)); %%% Labels are in second column (1st column is channel number)
    
    if size(depth_electrode_info, 2) > 2
        location = strjoin(depth_el_info{index, 3:end}, ' '); %%% The text file is space delimited, so join third to last column with spaces for loctation    
    end
    
end
    
end


%%% Function to determine whether subject is enrolled in UPenn experiments.
%%% Penn experiments contain .json files. So attempt to retrieve that data
%%% from these files and confirm that the entered ID is in R1999T format
%%% Would still need to confirm for every patient in case they didn't do
%%% Penn experiment sessions.

function [is_penn_subject, penn_ID] = determine_if_penn(subject, n_microelectrodes)
%%% In case it is not found
is_penn_subject = false(n_microelectrodes, 1);
penn_ID = repelem({'O'}, n_microelectrodes, 1);

%%% Search both raw_directory (Elemem recordings) and behavioral_directory (behavioral data) for .json that could have Penn subject data
raw_directory = dir(fullfile('/path/to/subject_files', subject, 'raw/**/experiment_config.json'));
behavioral_directory = dir(fullfile('/path/to/subject_files', subject, 'behavioral/**/session.json'));

if ~isempty(raw_directory)

    raw_directory = raw_directory(1);
    json_file = fullfile({raw_directory.folder}, {raw_directory.name});
    json_file = json_file{1};
    json_struct = jsondecode(fileread(json_file));
    subject = json_struct.subject;
    matched_subject_code = regexp(subject, 'R[0-9][0-9][0-9][0-9]T', 'match');
    
    if ~isempty(matched_subject_code)
        is_penn_subject = true(n_microelectrodes, 1);
        penn_ID = repelem(matched_subject_code(1), n_microelectrodes, 1);
    end
        
end

if ~isempty(behavioral_directory)

    behavioral_directory = behavioral_directory(1);
    json_file = fullfile({behavioral_directory.folder}, {behavioral_directory.name});
    json_file = json_file{1};
    json_text = fileread(json_file);
    matched_subject_code = regexp(json_text, 'R[0-9][0-9][0-9][0-9]T', 'match');
    
    if ~isempty(matched_subject_code)
        is_penn_subject = true(n_microelectrodes, 1);
        penn_ID = repelem(matched_subject_code(1), n_microelectrodes, 1);
    end
    
end

end