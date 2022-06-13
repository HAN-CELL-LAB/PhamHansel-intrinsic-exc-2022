clc; clear; close all; 
run start_up.m

%% Define data paths
main_data_path = 'data/gill-data/Layer II.III S1 BC IP'; 
DAT_path = fullfile(main_data_path, 'dat');
ITX_path = fullfile(main_data_path, 'itx'); 
INFO_path = fullfile(main_data_path, 'info');
MAT_path = fullfile(main_data_path, 'mat'); 
PROCESS_path = fullfile(main_data_path, 'processed'); 

%% File lists
ITX_files = dir(fullfile(ITX_path, '*.itx')); 
ITX_files = {ITX_files.name};

DAT_files = dir(fullfile(DAT_path, '*.dat')); 
DAT_files = {DAT_files.name};

ITX_and_DAT_files = [ITX_files, DAT_files]; 

INFO_file = fullfile(INFO_path, 'info-summary.xlsx');

%% Info sheet
info_table = readtable(INFO_file, 'Sheet', 'File');

num_cells = size(info_table, 1); 

data_file_patterns = struct;
data_file_patterns.DAT_files = '.dat';
data_file_patterns.ITX_files  = '.itx';
data_file_patterns.ITX_Base_files = 'Baseline.itx';
data_file_patterns.ITX_Post_files = 'Post.itx';
file_pat_files = fieldnames(data_file_patterns);

% insert more columns 
additional_columns = [{'full_ID'; 'have_data'}; file_pat_files; {'MAT_files'}]; 
for i = 1:length(additional_columns)
    col_name = additional_columns{i};
    info_table.(col_name) = arrayfun(@(x) '', ones(num_cells,1),'uni',0);
end
info_table.have_data = false(num_cells, 1);

%% Process files    

get_contains_of_list = @(file_list, pattern) file_list(contains(file_list, pattern, 'IgnoreCase', true));

t0 = tic; 
for i = 1:num_cells
    row_data = info_table(i,:);
    full_ID = [row_data.Prefix{1} row_data.ID{1}]; 
    
    row_data.full_ID = {full_ID};
    files_contain_prefix = ITX_and_DAT_files(contains(ITX_and_DAT_files, full_ID));

     files_based_on_types = structfun(@(pat) ...
        get_contains_of_list(files_contain_prefix,pat), ...
        data_file_patterns, 'uni', 0); 
    
    for j = 1:length(file_pat_files)
        pat = file_pat_files{j}; 
        file_name = files_based_on_types.(pat); 
        if isempty(file_name), continue; end
        if length(file_name) > 1, file_name = {file_name}; end
        row_data.(pat) = file_name; 
    end
    
    row_data.have_data = ~all(structfun(@isempty, files_based_on_types));
    
    row_ITX_files = row_data.ITX_files{1};
    if ~isempty(row_ITX_files) 
        itx2mat(row_ITX_files, ITX_path, MAT_path, 'N');
        row_data.MAT_files = {regexprep(row_ITX_files, '.itx$', '.mat')};
    end
    
    info_table(i,:) = row_data;
    
    t1 = toc(t0)/60; 
    fprintf('>>> Progress: %d/%d done. Elapsed %.1f min. ETA %.1f min.\n', ...
        i, num_cells, t1, num_cells*t1/i - t1);
end

%% Concatenate 
groupbase_table = cellfun(@(x) readtable(INFO_file, 'Sheet', x), {'Groups', 'Baseline'}, 'uni', 0);
groupbase_table = join(groupbase_table{:});

info_summary_table = join(info_table, groupbase_table);

%% Save 
save(fullfile(INFO_path, 'info-summary.mat'), 'info_summary_table', '*_path');