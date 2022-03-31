clc; clear; close all;
run start_up.m
run extract_gill_data_config.m

%% Define data paths
main_data_path = 'data/gill-data/Layer II.III S1 BC IP';
INFO_path = fullfile(main_data_path, 'info');

load(fullfile(INFO_path, 'info-summary.mat'));

%% Filter selection
% select_criteria = @(MAT_files, Group, Ihold) ...
%     all(contains(MAT_files{1}, {'Baseline', 'Post'}, 'IgnoreCase', true)) && ...
%     (strcmpi(Group, 'WT Somatic') || strcmpi(Group, 'WT Synaptic')) && ...
%     ~isnan(Ihold);
select_criteria = @(MAT_files, Group, Ihold) ...
    all(contains(MAT_files{1}, {'Baseline', 'Post'}, 'IgnoreCase', true)) && ...
    contains(Group, 'WT') && ...
    ~isnan(Ihold);

sel_col = rowfun(select_criteria,  info_summary_table, 'InputVariables', ...
    {'MAT_files', 'Group', 'Ihold'}, ...
    'OutputVariableNames', 'selected_for_process');

info_summary_table = [info_summary_table, sel_col];

%% Select for representative traces
represent_struct = struct; 
represent_struct.base = struct('portion', 1);
represent_struct.post = struct('portion', [0.2, 0.9]);

%% Analysis table
selected_info_table = info_summary_table(info_summary_table.selected_for_process,:);
num_cells = size(selected_info_table,1);

% dry run to get empty struct to initialize
init_info = table2struct(selected_info_table(1,:));
[init_sumanly, init_represent] = process_data_per_cell(init_info, MAT_path, configs, represent_struct, true);

analysis_table = table; 
analysis_table.full_ID = repmat({''}, [num_cells ,1]);
analysis_table.info = repmat(init_info, [num_cells ,1]);
analysis_table.analysis = repmat(init_sumanly, [num_cells ,1]);
analysis_table.representative = repmat(init_represent, [num_cells ,1]);
%%
t0 = tic;
for i = 1:num_cells
    cell_info = table2struct(selected_info_table(i,:));
    cell_fullID = cell_info.full_ID; 
    analysis_table(i,:).full_ID = {cell_fullID};
    analysis_table(i,:).info = cell_info;
    [summary_struct, represent_traces] = process_data_per_cell(cell_info, MAT_path, configs, represent_struct);
    analysis_table(i,:).analysis = summary_struct;
    analysis_table(i,:).representative = represent_traces;
    t1 = toc(t0)/60;
    fprintf('+ Cell "%s" done. Progress: %d/%d done. Elapsed %.1f min. ETA %.1f min.\n', ...
        cell_fullID, i, num_cells, t1, num_cells*t1/i - t1);
end

save(fullfile(PROCESS_path, 'analysis-summary.mat'), ...
    'analysis_table', '*path', 'info_summary_table');

%% Function to extract & process per cell/file
function [summary_struct, represent_traces] = process_data_per_cell(cell_info_struct, data_path, ...
    config_struct, represent_struct, dry_run)

if ~exist('dry_run', 'var')
    dry_run = false;
end

mat_files = fullfile(data_path, cell_info_struct.MAT_files);
base_file = mat_files{contains(mat_files, 'Base')};
post_file = mat_files{contains(mat_files, 'Post')};

info_struct = struct; 
info_struct.Ihold = cell_info_struct.Ihold; 

[base_summary, base_repr] = process_data_per_file(base_file, info_struct, config_struct, represent_struct.base, dry_run);
[post_summary, post_repr] = process_data_per_file(post_file, info_struct, config_struct, represent_struct.post, dry_run);

summary_time = [-length(base_summary):-1,1:length(post_summary)];

summary_struct = structarray_to_struct([base_summary; post_summary]);
summary_struct.time_vec = summary_time;

represent_traces = struct;
represent_traces.base = base_repr;
represent_traces.post = post_repr;

end

function [summary_struct, represent_traces] = process_data_per_file(...
    file_name, cell_info, config_struct, represent_struct, dry_run)
if ~exist('dry_run', 'var')
    dry_run = false;
end

load(file_name, 'recordings');

n_trace = length(recordings);
len_trace = length(recordings(1).data);
dt = recordings(1).dt * 1e3;
t = (0:len_trace-1)*dt;

represent_traces = struct('t', nan, 'Vm', nan);

init_src_struct = struct('configs', struct, 'spike_properties', struct,...
    'AP_properies', struct, 'passive_properties', struct); 

init_sum_struct = struct(...
    'num_spikes', nan, 'firing_rate', nan, 'latency', nan, ...
    'CV_ISI', nan, 'adapt_index', nan, ...
    'AP_Vm_mean', nan, 'fAHP_Vm_mean', nan, ...
    'Vthres_first', nan, 'Vthres_mean_afterfirst', nan, 'Vthres_mean_all', nan, ...
    'Vbase', nan, 'Vrest_1', nan, 'Vrest_2', nan, 'Rin', nan, ...
    'src', init_src_struct);

fields_from_spikeprops = {'num_spikes', 'firing_rate', 'latency', 'CV_ISI', 'adapt_index'}; 

summary_struct = repmat(init_sum_struct, [n_trace,1]);

if dry_run
    return;
end

spike_params = config_struct.spike_params;
stim_info = config_struct.stim_info;
AP_params = config_struct.AP_params;
pass_params = config_struct.pass_params;


type_of_reprselect = fieldnames(represent_struct);
if length(type_of_reprselect) ~= 1
    error('Allow only one field for representative trace struct');
end
type_of_reprselect = type_of_reprselect{1}; 
ind_savedtraces = represent_struct.(type_of_reprselect); 

switch upper(type_of_reprselect)
    case {'NUM', 'IND'}
        ind_savedtraces = round(ind_savedtraces);
    case {'PORTION','P','PERCENT'}
        ind_savedtraces = round(ind_savedtraces*n_trace);
    otherwise
        error('"%s" is not an accepted type of select for representative traces', type_of_reprselect);
end
ind_savedtraces = bound_minmax(ind_savedtraces, 1, n_trace); 
n_savedtraces = length(ind_savedtraces);
represent_traces.t = t; 
represent_traces.Vm = zeros(len_trace, n_savedtraces);
cnt_savedtraces = 1; 

for i = 1:n_trace
    trace_struct = recordings(i);
    note_msg = '';
    if isfield(trace_struct, 'note_msg')
        note_msg = trace_struct.note_msg;
    end
    
    if ~isempty(note_msg)
        continue;
    end
    
    Vm = trace_struct.data * 1e3;
    
    spike_properties = calculate_spike_properties(t, Vm, spike_params, stim_info);
    
    spike_indices = spike_properties.spike_inds;
    if ~isempty(spike_indices)
        AP_properies = arrayfun(@(AP_ind) calculate_AP_properties(Vm, AP_ind, dt, AP_params), ...
            spike_indices, 'uni', 1);
        AP_properies = structarray_to_struct(AP_properies);
    else
        AP_properies =  calculate_AP_properties([],[],[],[],true); % dry run
    end
    
    passive_properties = calculate_passive_properties(t, Vm, cell_info, stim_info, pass_params);
    
    % start saving     
    for j = 1:length(fields_from_spikeprops)
        field_j = fields_from_spikeprops{j}; 
        summary_struct(i).(field_j) = spike_properties.(field_j);
    end
    
    summary_struct(i).Vbase = passive_properties.Vbase;
    summary_struct(i).Rin = passive_properties.Rin;
    summary_struct(i).Vrest_1 = passive_properties.Vrest_1;
    summary_struct(i).Vrest_2 = passive_properties.Vrest_2;
    
    summary_struct(i).AP_Vm_mean = mean(AP_properies.AP_Vm);
    summary_struct(i).fAHP_Vm_mean = mean(AP_properies.fAHP_Vm);
    
    Vthres = AP_properies.AP_thres_Vm;
    if ~isempty(Vthres) && all(~isnan(Vthres))
        summary_struct(i).Vthres_first = Vthres(1);
        summary_struct(i).Vthres_mean_all = mean(Vthres);
    end
    if length(Vthres) > 1
        summary_struct(i).Vthres_mean_afterfirst = mean(Vthres(2:end));
    end
    
    summary_struct(i).src.configs = config_struct;
    summary_struct(i).src.spike_properties = spike_properties;
    summary_struct(i).src.AP_properies = AP_properies;
    summary_struct(i).src.passive_properties = passive_properties;      
    
    if any(ind_savedtraces == i)
        represent_traces.Vm(:,cnt_savedtraces) = Vm; 
        cnt_savedtraces = cnt_savedtraces + 1;
    end
end

end