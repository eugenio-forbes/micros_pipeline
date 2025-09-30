%%% Base function of Micros Pipeline, Step 7: Manual Review
%%% This function will generate job lists that can be used
%%% in Combinato's GUI to perform a review of the spike
%%% sorting results, which can be modified in the GUI.
%%% This speeds up the process of review as it eliminated the
%%% need of individually loading files for review.

function n07_make_manual_list(varargin)
if isempty(varargin)   %%% Will make a list of all combinato directories
    %%% Specify where to save do_manual_neg/pos.txt file%%%
    
    root_directory = '/path/to/micros_pipeline/parent_directory';
    review_type = 'all'; %%% See switch-case statement for more information on types.
    
    data_directory = fullfile(root_directory, 'micros_database');
    save_directory = fullfile(root_directory, 'micros_pipeline/process_files/n07_manual');
    
    load(fullfile(data_directory, 'progress_table.mat'), 'progress_table');
    clustered_recordings = progress_table(progress_table.clustered, :); clear progress_table;
    
    switch(review_type)
        case 'all'
            save_directory = fullfile(save_directory, 'review_lists');  %%% Divided by task. All micros listed.
        case 'subject_lists'
            save_directory = fullfile(save_directory, 'subject_lists'); %%% Divided by subject. All micros listed.
        case 'manual_lists'
            save_directory = fullfile(save_directory, 'manual_lists');  %%% Divided by task. Only micros missing manual review listed.
    end
    
else   %%% Will make a list of combinato directories depending on type (varargin(2)). If selective for manual review process list, takes vargin(3) for list of micros.
    root_directory = varargin{1};
    review_type = varargin{2};
    
    data_directory = fullfile(root_directory, 'micros_database');
    save_directory = fullfile(root_directory, 'micros_pipeline/process_files/n07_manual');
    
    switch(review_type)
        case 'all'
            load(fullfile(data_directory, 'progress_table.mat'), 'progress_table');
            clustered_recordings = progress_table(progress_table.clustered, :); clear progress_table;
            save_directory = fullfile(save_directory, 'review_lists');  %%% Divided by task. All micros listed.
        case 'subject_lists'
            load(fullfile(data_directory, 'progress_table.mat'), 'progress_table');
            clustered_recordings = progress_table(progress_table.clustered, :); clear progress_table;
            save_directory = fullfile(save_directory, 'subject_lists'); %%% Divided by subject. All micros listed.
        case 'manual_lists'
            save_directory = fullfile(save_directory, 'manual_lists');  %%% Divided by task. Only micros missing manual review listed.
            clustered_recordings = varargin{3};
    end
end

%%% List all tasks with behavioral events
tasks = {'all', 'AR', 'AR-scopolamine', 'AR-stim', 'BN-stim', 'Courier', 'FR', 'FR-scopolamine', 'no-task', 'SR'};

if ~contains(review_type, 'subject')
    for idx = 1:length(tasks)
        paths = {};
        
        task_folder = tasks{idx};
        
        switch task_folder
            case 'all'
                these_recordings = clustered_recordings;
            case 'BN-stim'
                these_recordings = clustered_recordings(contains(clustered_recordings.task, task_folder), :);
            otherwise
                these_recordings = clustered_recordings(strcmp(clustered_recordings.task, task_folder), :);
        end
        
        for mdx = 1:height(these_recordings)
        
            subject = these_recordings(mdx, :).subject{:};
            date = these_recordings(mdx, :).date{:};
            task = these_recordings(mdx, :).task{:};
            part = these_recordings(mdx, :).part;
            folder = strcat(date, '_', task);
            if part > 0
                folder = strcat(folder, sprintf('_part%d', part));
            end
            bank = these_recordings(mdx, :).bank{:};
            
            combinato_folder = fullfile(data_directory, subject, folder, 'combinato_files', sprintf('Bank%s', bank));
            
            sort_paths = dir(fullfile(combinato_folder, 'NS6*', 'data_NS6*.h5'));
            sort_paths = fullfile({sort_paths.folder}, {sort_paths.name});
            
            paths = [paths; sort_paths'];
        end
        
        this_save_directory = fullfile(save_directory, task_folder);
        manual_file1 = fullfile(this_save_directory, 'do_manual_neg.txt');
        manual_file2 = fullfile(this_save_directory, 'do_manual_pos.txt');
        
        if isempty(these_recordings)
            delete(manual_file1);
            delete(manual_file2);
        else
            if ~exist(this_save_directory, 'dir')
                mkdir(this_save_directory);
            end
            
            if isfile(manual_file1)
                delete(manual_file1);
            end
            if isfile(manual_file2)
                delete(manual_file2);
            end
            
            writecell(paths, manual_file1);
            writecell(paths, manual_file2);
        end
    end
else
    subjects = unique(clustered_recordings.subject);
    
    for idx = 1:length(subjects)
    
        paths = {};
        
        subject = subjects{idx};
        these_recordings = clustered_recordings(contains(clustered_recordings.subject, subject), :);
        
        for mdx = 1:height(these_recordings)
        
            date = these_recordings(mdx, :).date{:};
            task = these_recordings(mdx, :).task{:};
            part = these_recordings(mdx, :).part;
            folder = strcat(date, '_', task);
            if part > 0
                folder = strcat(folder, sprintf('_part%d', part));
            end
            bank = these_recordings(mdx, :).bank{:};
            
            combinato_folder = fullfile(data_directory, subject, folder, 'combinato_files', sprintf('Bank%s', bank));
            
            sort_paths = dir(fullfile(combinato_folder, 'NS6*', 'data_NS6*.h5'));
            sort_paths = fullfile({sort_paths.folder}, {sort_paths.name});
            
            paths = [paths; sort_paths'];
        end
        
        this_save_directory = fullfile(save_directory, subject);
        manual_file1 = fullfile(this_save_directory, 'do_manual_neg.txt');
        manual_file2 = fullfile(this_save_directory, 'do_manual_pos.txt');
        
        if isempty(these_recordings)
            delete(manual_file1);
            delete(manual_file2);
        else
        
            if ~exist(this_save_directory, 'dir')
                mkdir(this_save_directory);
            end
            
            if isfile(manual_file1)
                delete(manual_file1);
            end
            if isfile(manual_file2)
                delete(manual_file2);
            end
            
            writecell(paths, manual_file1);
            writecell(paths, manual_file2);
        end
    end  
end
end
