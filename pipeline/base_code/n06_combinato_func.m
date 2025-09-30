%%% Base function of Micros Pipeline, Step 6: Clustering
%%% This function will execute python code using system
%%% shells to cluster spike data using Combinato.
%%% It may use 'standard' parameters for clustering
%%% or loop through a combination of parameters to perform
%%% grid search of ideal parameters for clustering.

function [error_flag, error_message] = n06_combinato_func(varargin)
if isempty(varargin)                                               %%% To run manually edit values below
    root_directory = '/path/to/micros_pipeline/parent_directory';  %%% Root directory with pipeline and database folders
    subject = 'SC000';                                             %%% Subject code in SC000 format
    folder = 'yyyy-mm-dd_task-code_part1';                         %%% Folder in yyyy-mm-dd_task-code format, or yyyy-mm-dd_task-code_part1 format if more than 1 part. 
    bank = 'A';                                                    %%% Recording hardware bank character ('A', 'B', 'C', 'D')
    clustering_mode = 'standard';                                  %%% Clustering modes:
                                                                   %%% 'standard' uses custom_options.txt only, 
                                                                   %%% 'grid-search' uses combinations of parameters from template_options.txt in addition to custom_options.txt
else                                                               %%% Otherwise this is the order they should be entered into function, following above format
    root_directory = varargin{1};
    subject = varargin{2};
    folder = varargin{3};
    bank = varargin{4};
    clustering_mode = varargin{5};
end

%%% List directories
data_directory = fullfile(root_directory, 'micros_database', subject, folder, 'combinato_files', sprintf('Bank%s', bank));
delete_previous_attempts(data_directory);    %%% Because duplicate jobs cause errors
template_directory = fullfile(root_directory, 'micros_pipeline/pipeline/base_code/combinato_templates');

options_text_file =  fullfile(template_directory, 'custom_options.txt');
options_text = fileread(options_text_file);

options_template_file =  fullfile(template_directory, 'template_options.txt');
options_template = fileread(options_template_file);

extraction_template_file = fullfile(template_directory, 'template_extract.txt');
extraction_template = fileread(extraction_template_file);

cluster_template_file = fullfile(template_directory, 'template_cluster.txt');
cluster_template = fileread(cluster_template_file);

%%% Print custom options to local_options.py file and make executable
options_file = fullfile(data_directory, 'local_options.py');
file_id = fopen(options_file, 'w');
fprintf(file_id, options_text);
fclose(file_id);
system(sprintf('chmod +x %s', options_file));

%%% Print bash file for extracting spikes and make executable
extraction_text = sprintf(extraction_template, data_directory); %%% Modifying clustering directory to a bank of a session
extraction_file = fullfile(data_directory, 'extract.sh');
file_id = fopen(extraction_file, 'w');
fprintf(file_id, extraction_text);
fclose(file_id);
system(sprintf('chmod +x %s', extraction_file));

%%% Values to be returned by function
error_flag = 0;
error_message = '';

%%% Run extraction of spikes
[~, command_output] = system(sprintf('bash %s', extraction_file));
            
%%% End function execution if combinato outputs traceback or error
%%% Return combinato output and flag error
if contains(command_output, {'traceback', 'error'}, 'IgnoreCase', true)
    error_flag = 1;
    error_message = command_output;
    return
else
    delete(extraction_file)
end

%%% Entries for clustering template
entries = [];
entries = [entries; {data_directory}];       %%% data_directory first
label = 'options00';
entries = [entries; repelem({label}, 6, 1)]; %%% then 6 times label

%%% Print bash file for clustering spikes and make executable
cluster_text = sprintf(cluster_template, entries{:}); %%% Input entries to template
cluster_file = fullfile(data_directory, 'cluster.sh');
file_id = fopen(cluster_file, 'w');
fprintf(file_id, cluster_text); fclose(file_id);
system(sprintf('chmod +x %s', cluster_file));

%%% Run clustering in system shell
[~, command_output] = system(sprintf('bash %s', cluster_file));

%%% End function execution if combinato outputs traceback or error
%%% Return combinato output and flag error
if contains(command_output, {'traceback', 'error'}, 'IgnoreCase', true)
    error_flag = 1;
    error_message = command_output;
    return
else
    if isfile(cluster_file)
        delete(cluster_file);
    end
    if isfile(cluster_file)
        delete(options_file);
    end
    pycache_directory = fullfile(data_directory, '__pycache__');
    if isfolder(pycache_directory)
        rmdir(pycache_directory, 's');
    end
end

if strcmp(clustering_mode, 'grid-search')
    %%% Get combinations of parameters to sort
    grid_parameters = get_grid_parameters();
    
    %%% Loop through parameters and run sorting
    for idx = 1:height(grid_parameters)
        
        %%% Entries for option template
        this_parameters = grid_parameters{idx, :};
        entries = [this_parameters{:}];
        
        %%% Print custom options to local_options.py file and make executable
        options_text = sprintf(options_template, entries(:));
        options_file = fullfile(data_directory, 'local_options.py');
        file_id = fopen(options_file, 'w'); 
        fprintf(file_id, options_text); 
        fclose(file_id);
        system(sprintf('chmod +x %s', options_file));
        
        %%% Entries for clustering template
        entries = [];
        entries = [entries; {data_directory}];       %%% data_directory first
        label = sprintf('options%02d', idx);
        entries = [entries; repelem({label}, 6, 1)]; %%% then 6 times label
        
        %%% Print bash file for clustering spikes and make executable
        cluster_text = sprintf(cluster_template, entries{:}); %%% Input entries to template
        cluster_file = fullfile(data_directory, 'cluster.sh');
        file_id = fopen(cluster_file, 'w');
        fprintf(file_id, cluster_text);
        fclose(file_id);
        system(sprintf('chmod +x %s', cluster_file));
        
        %%% Run clustering in system shell
        [~, command_output] = system(sprintf('bash %s', cluster_file));
        
        %%% End function execution if combinato outputs traceback or error
        %%% Return combinato output and flag error
        if contains(command_output, {'traceback', 'error'}, 'IgnoreCase', true)
            error_flag = 1;
            error_message = command_output;
            return
        else
            delete(cluster_file)
            delete(options_file)
            rmdir(fullfile(data_directory, '__pycache__'), 's');
        end
    end
end

end


function grid_parameters = get_grid_parameters()
%%% Define sets of parameters to compare

artifact = [{[4, 1.8, 1.3, 0.9]}; {[5, 2, 1.5, 0.9]}];
%%%   [number of maxima, maxima ratio, max min ratio, peak to peak ratio] 
%%%   [{strict}; vs {default}]

groups = [{1.65}; {1.95}];
%%%   'MaxDistMatchGrouping'
%%%   [{tight}; vs {loose}]

clusters = [{[15, 0.65, 2.6, 8]}; {[20, 0.85, 3.4, 5]}];
%%%   [MinSpikesPerClusterMultiSelect, FirstMatchFactor, SecondMatchFactor, MaxClustersPerTemp]
%%%   [{tight}; vs {loose}]

recursion = [{2}; {1}]; %%% Recluster big clusters vs not

%%% Make table with all possible combinations of sets. Should be 2^4.
[Ax, Bx, Cx, Dx] = ndgrid(1:numel(artifact), 1:numel(groups), 1:numel(clusters), 1:numel(recursion));        
grid_parameters = table;
grid_parameters.artifact = artifact(Ax(:));
grid_parameters.groups = groups(Bx(:));
grid_parameters.clusters = clusters(Cx(:));
grid_parameters.recursion = recursion(Dx(:));

end


function delete_previous_attempts(data_directory)
    %%% Delete previous attempt files if there are any. Repeat run would mess things up.
    
    invalid = {'.', '..'};
    files_to_keep = {'.hdf5'};
    
    directory_list = dir(data_directory);
    directory_list = directory_list(~ismember({directory_list.name}, invalid));
    
    files = directory_list(~[directory_list.isdir]);
    files = files(~contains({files.name}, files_to_keep));
    combinato_files = fullfile({files.folder}, {files.name});
    
    for idx = 1:length(combinato_files)
        this_file = combinato_files{idx};
        delete(this_file);
    end
    
    directories = directory_list([directory_list.isdir]);
    combinato_folders = fullfile({directories.folder}, {directories.name});
    for idx = 1:length(combinato_folders)
        this_folder = combinato_folders{idx};
        rmdir(this_folder, 's');
    end
end
