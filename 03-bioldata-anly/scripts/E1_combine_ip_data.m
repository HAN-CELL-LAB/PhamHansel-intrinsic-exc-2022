clc; clear; close all;
run start_up.m

data_sources = {...
    'data/gill-data/Layer II.III S1 BC IP/processed/analysis-summary.mat',...
    'data/supp-data/processed/analysis-summary.mat'};

selected_info_fields = {'Prefix', 'full_ID', 'ID', 'Group'};
selected_anly_fields = {'num_spikes', 'Rin', 'Vbase', 'Vrest_1', 'Vrest_2', ...
                        'Vthres_first', 'Vthres_mean_afterfirst', 'Vthres_mean_all', ...
                        'time_vec'};
                    
filter_fn = @(T) filter_struct_inside_table(T, ...
    'info', selected_info_fields, ...
    'analysis', selected_anly_fields);

data_tables = cellfun(@(f) load_and_filter_data(f, filter_fn), data_sources, 'uni', 0);
analysis_table = vertcat(data_tables{:});
analysis_table = rename_table_groups(analysis_table);

save('data/ip-data/analysis-summary.mat', 'analysis_table', 'data_sources');

function T = rename_table_groups(T)
    for i = 1:height(T)
        T(i,:).info.Group = rename_exp_group(T(i,:).info.Group);
    end
end

function s = rename_exp_group(s)
    % very specific to this 
    s = strip(lower(strrep(s, 'WT', '')));
    is_elec = contains(s, 'somatic') || contains(s, 'synaptic');
    is_chol = contains(s, 'oxo-m') || contains(s, 'oxom');

    if is_chol && is_elec
        s = 'cholinergic paired';
        return;
    end

    if is_elec
        s = 'electric';
    end

    if is_chol
        s = 'cholinergic';
    end
    
    
end

function T = load_and_filter_data(file_name, filter_fn)
    T = load(file_name).analysis_table; 
    T = filter_fn(T); 
end


