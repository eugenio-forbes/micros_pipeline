function reset_session(varargin)
if isempty(varargin) %%% To run from editor edit these values
    root_directory = '/project/TIBIR/Lega_lab/shared/lega_ansir/Single_Units_BR';
    target_subject = 'SC000';
    target_date = 'yyyy-mm-dd';
    target_task = 'task-name';
    target_part = 0;
else
    root_directory = varargin{1};
    target_subject = varargin{2};
    target_date = varargin{3};
    target_task = varargin{4};
    target_part = varargin{5};
end

data_directory = fullfile(root_directory, 'micros_database');
progress_table_file = fullfile(data_directory, 'progress_table.mat');
load(progress_table_file, 'progress_table');

has_subject = strcmp(progress_table.subject, target_subject);
has_date = strcmp(progress_table.date, target_date);
has_task = strcmp(progress_table.task, target_task);
has_part = progress_table.part == target_part;

indices = find(has_subject & has_date & has_task & has_part);

folder = strcat(target_date, '_', target_task);
if target_part > 0
    folder = strcat(folder, '_', sprintf('part%d', target_part));
end

if ~isempty(indices)
    
    progress_table.split(indices) = false(length(indices), 1);
    progress_table.rerefed(indices) = false(length(indices), 1);
    progress_table.clean(indices) = false(length(indices), 1);
    progress_table.rescaled(indices) = false(length(indices), 1);
    progress_table.hdf5(indices) = false(length(indices), 1);
    progress_table.clustered(indices) = false(length(indices), 1);
    progress_table.manual(indices) = false(length(indices), 1);
    progress_table.spikes_timed(indices) = false(length(indices), 1);
    progress_table.modal(indices) = false(length(indices), 1);
    progress_table.events(indices) = false(length(indices), 1);
    progress_table.has_events(indices) = false(length(indices), 1);
    progress_table.has_sync(indices) = false(length(indices), 1);
    progress_table.aligned(indices) = false(length(indices), 1);
    
    save(progress_table_file, 'progress_table');
    
    this_directory = fullfile(data_directory, target_subject, folder);
    split_directory = fullfile(this_directory, 'split');
    reref_directory = fullfile(this_directory, 'reref');
    clean_directory = fullfile(this_directory, 'clean');
    rescale_directory = fullfile(this_directory, 'rescaled');
    combinato_directory = fullfile(this_directory, 'combinato_files');
    lfp_directory = fullfile(this_directory, 'lfp');
    modal_directory = fullfile(this_directory, 'modal');
    align_directory = fullfile(this_directory, 'alignment_report');
    
    if isfolder(split_directory)
        rmdir(split_directory, 's');
    end
    
    if isfolder(reref_directory)
        rmdir(reref_directory, 's');
    end
    
    if isfolder(clean_directory)
        rmdir(clean_directory, 's');
    end
    
    if isfolder(rescale_directory)
        rmdir(rescale_directory, 's');
    end
    
    if isfolder(combinato_directory)
        rmdir(combinato_directory, 's');
    end
    
    if isfolder(lfp_directory)
        rmdir(lfp_directory, 's');
    end
    
    if isfolder(modal_directory)
        rmdir(modal_directory, 's');
    end
    
    if isfolder(align_directory)
        rmdir(align_directory, 's');
    end

    reref_plots = sprintf(fullfile(root_directory, 'micros_pipeline/process_files/n02_reref/plots/%s/%s*'), target_subject, folder);
    delete(reref_plots)
    
    clean_plots = sprintf(fullfile(root_directory, 'micros_pipeline/process_files/n03_clean/plots/%s'), target_subject, folder);
    delete(clean_plots)
    
    error_logs = sprintf(fullfile(root_directory, 'micros_pipeline/error_logs/%s/%s'), target_subject, folder);
    if isfolder(error_logs)
        rmdir(error_logs, 's');
    end
    
    fprintf('Target session %s from subject %s has been reset.\n', folder, target_subject)    
else
    fprintf('Target session %s from subject %s not logged or non-existent. Check table and try again.\n', folder, target_subject)    
end

end