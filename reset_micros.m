function reset_micros(varargin)
if isempty(varargin) %To run from editor edit these values
    root_directory = '/project/TIBIR/Lega_lab/shared/lega_ansir/Single_Units_BR';
else
    root_directory = varargin{1};
end

data_directory = fullfile(root_directory, 'micros_database');
progress_table_file = fullfile(data_directory, 'progress_table.mat');
load(progress_table_file, 'progress_table');

progress_table.split = false(height(progress_table), 1);
progress_table.rerefed = false(height(progress_table), 1);
progress_table.clean = false(height(progress_table), 1);
progress_table.rescaled = false(height(progress_table), 1);
progress_table.hdf5 = false(height(progress_table), 1);
progress_table.clustered = false(height(progress_table), 1);
progress_table.manual = false(height(progress_table), 1);
progress_table.spikes_timed = false(height(progress_table), 1);
progress_table.lfp = false(height_porgress_table, 1);
progress_table.modal = false(height(progress_table), 1);
progress_table.events = false(height(progress_table), 1);
progress_table.has_events = false(height(progress_table), 1);
progress_table.has_sync = false(height(progress_table), 1);
progress_table.aligned = false(height(progress_table), 1);

save(progress_table_file, 'progress_table');

for idx = 1:11
    if idx ~= 7
        delete(fullfile(root_directory, sprintf('micros_pipeline/process_files/n%02d_*/locks/*', idx)));
    end        
    delete(fullfile(root_directory, sprintf('micros_pipeline/process_files/n%02d_*/check/*', idx)));
end

reref_plots = fullfile(root_directory, 'micros_pipeline/process_files/n02_reref/plots');
if isfolder(reref_plots)
    rmdir(reref_plots, 's');
end

clean_plots = fullfile(root_directory, 'micros_pipeline/process_files/n03_clean/plots');
if isfolder(clean_plots)
    rmdir(clean_plots, 's');
end

error_logs = fullfile(root_directory, 'micros_pipeline/error_logs');
if isfolder(error_logs)
    rmdir(error_logs, 's');
end

[~, unique_idx, ~] = unique(fullfile(progress_table.subject, progress_table.date, progress_table.task, num2str(progress_table.part)));

progress_table = progress_table(unique_idx, :);

for idx = 1:height(progress_table)
    
    subject = progress_table(idx, :).subject{:};
    date = progress_table(idx, :).date{:};
    task = progress_table(idx, :).task{:};
    part = progress_table(idx, :).part;
    folder = strcat(date, '_', task);
    
    if part > 0
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

end