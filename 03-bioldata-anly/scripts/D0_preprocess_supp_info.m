clc; clear; close all; 
run start_up.m

%% Define data paths
main_data_path = 'data/supp-data/'; 
DAT_path = fullfile(main_data_path, 'dat');
ITX_path = fullfile(main_data_path, 'itx'); 
INFO_path = fullfile(main_data_path, 'info');
MAT_path = fullfile(main_data_path, 'mat'); 
PROCESS_path = fullfile(main_data_path, 'processed'); 

%% Info file
INFO_xls_file = fullfile(INFO_path, 'info-summary.xlsx');
INFO_mat_file = fullfile(INFO_path, 'info-summary.mat');

info_summary_table_xls = load_info(INFO_xls_file);
info_summary_table_mat = load(INFO_mat_file).info_summary_table;
proc_rows = find(info_summary_table_mat.processed);
if isempty(proc_rows)
    info_summary_table = info_summary_table_xls;
else
    proc_old_rows = info_summary_table_mat(proc_rows, :);
    proc_old_IDs = proc_old_rows.ID;
    
    new_rows = cellfun(@(x) all(~contains(proc_old_IDs, x)), info_summary_table_xls.ID);
    new_rows = info_summary_table_xls(new_rows, :);
    
    info_summary_table = [proc_old_rows; new_rows];
end

%% Save Info
save(INFO_mat_file, 'info_summary_table', '*_path');

%% Extract specific base and post protocols 
unprocessed_rows = info_summary_table(~info_summary_table.processed,:);

for i = 1:height(unprocessed_rows)
    row_info = table2struct(unprocessed_rows(i,:));
    save_base_and_post(MAT_path, row_info);    
end


%% Helper
function tbl = load_info(file_name)
    special_cols = {'Base_ID','Induction_ID','Post_ID', 'I_hold_base', 'I_hold_post'};
    opts = detectImportOptions(file_name);
    opts = setvartype(opts, special_cols, 'string');
    tbl = readtable(file_name, opts);

    for i = 1:length(special_cols)
        col = special_cols{i};
        tbl.(col) = cellfun(@str2num, tbl.(col), 'uni', false);
    end

end


function save_base_and_post(MAT_path, info)
    mat_file_name = fullfile(MAT_path, info.MAT_file); 
    base_id_vec = info.Base_ID; 
    post_id_vec = info.Post_ID;
    I_hold_base = info.I_hold_base;
    I_hold_post = info.I_hold_post; 
    
    if length(base_id_vec) ~= length(I_hold_base) || ...
        length(post_id_vec) ~= length(I_hold_post) 
        error('The base/post ID vec and base/post I_hold vec need to have same length');
    end 
    
    tbl = load(mat_file_name).data_table;

    prot_indices = [base_id_vec, post_id_vec];
    prot_tags = [repmat({'base'},[1,length(base_id_vec)]), repmat({'post'}, [1,length(post_id_vec)])];
    I_holds = [I_hold_base, I_hold_post]; 
%     I_holds = ones(1, length(prot_indices)) * info.I_hold; 
    
    recordings = load_data(tbl, prot_tags, prot_indices, I_holds);
    
    save_file_name = fullfile(MAT_path, [info.full_ID '-basepost.mat']);
    save(save_file_name, 'recordings', 'info');
    
    fprintf('Cell "%s" baseline and post extracted at "%s"\n', info.full_ID, save_file_name);
end

function recordings = load_data(tbl, tags, indices, I_holds, channel)
    if ~exist('channel', 'var')
        channel = 'Vmon';
    end
    
    if length(tags) ~= length(indices)
        error('The "tags" and "indices" vector need to be of same length');
    end

    num_prots = length(tags);
    recordings = cell(1, length(indices));

    for i = 1:num_prots
        dat = table2struct(tbl(indices(i),:)); 
        Fs = dat.SR;
        dt = 1/Fs;

        sel_ind = cellfun(@(x) startsWith(x, channel), dat.ChName);
        if sum(sel_ind) ~= 1
            error('The selected channel data should only be 1 per protocol');
        end

        data_vec = dat.dataRaw{sel_ind}; 
        recordings{i} = struct(...
            'Fs', Fs, ...
            'dt', dt, ...
            'tag', tags{i}, ...
            'Ihold', I_holds(i), ...
            'data', mat2cell_splitdim(data_vec,2)' ...
        ); 
    end

    recordings = vertcat(recordings{:});

end