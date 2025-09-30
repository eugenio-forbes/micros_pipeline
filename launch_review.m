function launch_review(varargin)
if isempty(varargin)
    root_directory = '/path/to/micros_pipeline/parent_directory';
    review_type = 'manual';
    folder = 'all';
else
    root_directory = varargin{1};
    review_type = varargin{2};
    folder = varargin{3};
end

%%% List directories
manual_directory = fullfile(root_directory, 'micros_pipeline/process_files/n07_manual');
combinato_directory = fullfile(manual_directory, 'combinato_files');
switch review_type
    case 'manual'
        list_directory = fullfile(manual_directory, 'manual_lists');
    case 'task'
        list_directory = fullfile(manual_directory, 'review_lists');
    case 'subject'
        list_directory = fullfile(manual_directory, 'subject_lists');
end

available_folders = dir(list_directory);
available_folders = available_folders([available_folders.isdir]);
available_folders = {available_folders.name};
available_folders = available_folders(~ismember(available_folders, {'.', '..'}));

neg_file = 'do_manual_neg.txt';
pos_file = 'do_manual_pos.txt';

if ismember(folder, available_folders)
    list_directory = fullfile(list_directory, folder);
    if isfile(fullfile(combinato_directory, neg_file))
        delete(fullfile(combinato_directory, neg_file));
    end
    if isfile(fullfile(combinato_directory, pos_file))
        delete(fullfile(combinato_directory, pos_file));
    end
    if isfile(fullfile(list_directory, neg_file))
        copyfile(fullfile(list_directory, neg_file), fullfile(combinato_directory, neg_file));
    else
        fprintf('Jobs list for negative spikes in %s lists for %s folder is not available. Could not change jobs list.\n', review_type, folder);
    end
    if isfile(fullfile(list_directory, pos_file))
        copyfile(fullfile(list_directory, pos_file), fullfile(combinato_directory, pos_file));
    else
        fprintf('Jobs list for positive spikes in %s lists for %s folder is not available. Could not change jobs list.\n', review_type, folder);
    end
else
    fprintf('%s folder for %s is not available. Could not change combinato jobs lists.\n', review_type, folder);
end

end