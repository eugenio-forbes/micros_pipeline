%%% Base function of Micros Pipeline, Step 5: Converting to HDF5
%%% This function saves a copy of rescaled signal in HDF5
%%% format for every channel in a microelectrode.

function n05_hdf5_func(varargin)
if isempty(varargin)                                               %%% To run manually edit values below
    root_directory = '/path/to/micros_pipeline/parent_directory';  %%% Root directory with pipeline and database folders
    subject = 'SC000';                                             %%% Subject code in SC000 format
    folder = 'yyyy-mm-dd_task-code_part1';                         %%% Folder in yyyy-mm-dd_task-code format, or yyyy-mm-dd_task-code_part1 format if more than 1 part. 
    bank = 'A';                                                    %%% Recording hardware bank character ('A', 'B', 'C', 'D')
else                                                               %%% Otherwise this is the order they should be entered into function, following above format
    root_directory = varargin{1};
    subject = varargin{2};
    folder = varargin{3};
    bank = varargin{4};
end

%%% List directories
data_directory = fullfile(root_directory, 'micros_database', subject, folder, 'rescaled', sprintf('Bank%s', bank));
save_directory = strrep(data_directory, 'rescaled', 'combinato_files');
if ~exist(save_directory, 'dir')
    mkdir(save_directory);
end

file_paths = dir(fullfile(data_directory, '*.mat'));
file_names = {file_paths.name};
file_names = file_names(~contains(file_names, 'noise_stds'));

for idx = 1:length(file_names)
    
    this_file = file_names{idx};
    
    this_hdf5 = fullfile(save_directory, strrep(this_file, '.mat', '.hdf5'));
    
    this_file = fullfile(data_directory, this_file);
    
    data_struct = load(this_file);
    data_cell = struct2cell(data_struct);
    
    field_names = fieldnames(data_struct); clear data_struct
    
    if isfile(this_hdf5)
        delete(this_hdf5);
    end
    
    for idx = 1:numel(field_names)
        field_data = data_cell{idx};
        data_size = length(field_data);
        
        h5create(this_hdf5, ['/', field_names{idx}], data_size, 'Datatype', class(field_data));
        h5write(this_hdf5, ['/', field_names{idx}], field_data);
    end
    
    clear data_struct data_cell
end

end