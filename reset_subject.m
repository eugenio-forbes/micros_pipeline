function reset_subject(varargin)
if isempty(varargin) %%% To run from editor edit these values
    root_directory = '/path/to/micros_pipeline/parent_directory';
    target_subject = 'SC000';
else
    root_directory = varargin{1};
    target_subject = varargin{2};
end

data_directory = fullfile(root_directory, 'micros_database');
progress_table_file = fullfile(data_directory, 'progress_table.mat');
load(progress_table_file, 'progress_table');

indices = find(strcmp(progress_table.subject, target_subject));

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
    progress_table.has_sync_pulses(indices) = false(length(indices), 1);
    progress_table.aligned(indices) = false(length(indices), 1);
    
    save(progress_table_file, 'progress_table');
    
    subject_recordings = progress_table(indices, :);

    for idx = 1:height(subject_recordings)
        
        subject = subject_recordings(idx, :).subject{:};
        date = subject_recordings(idx, :).date{:};
        task = subject_recordings(idx, :).task{:};
        part = subject_recordings(idx, :).part;
        folder = strcat(date, '_', task);
        if part>0
            folder = strcat(folder, '_', sprintf('part%d', part));
        end
        
        this_directory = fullfile(data_directory, subject, folder);
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
    end
    
    reref_plots = sprintf(fullfile(root_directory, 'micros_pipeline/process_files/n02_reref/plots/%s'), target_subject);
    if isfolder(reref_plots)
        rmdir(reref_plots, 's');
    end
    
    clean_plots = sprintf(fullfile(root_directory, 'micros_pipeline/process_files/n03_clean/plots/%s'), target_subject);
    if isfolder(clean_plots)
        rmdir(clean_plots, 's');
    end
    
    error_logs = sprintf(fullfile(root_directory, 'micros_pipeline/error_logs/%s'), target_subject);
    if isfolder(error_logs)
        rmdir(error_logs, 's');
    end
    
    fprintf('Target subject %s has been reset.\n', target_subject)    
else
    fprintf('Target subject %s not logged or non-existent. Check table and try again.\n', target_subject)    
end

end