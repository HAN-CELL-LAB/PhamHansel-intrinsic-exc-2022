clc; clear; close all;
run start_up.m

%% Define data paths
main_data_path = 'data/pham-data/';
debug_fig_path = 'figures/pham-data/debug';

mat_path = fullfile(main_data_path, 'mat');
proc_path = fullfile(main_data_path, 'processed'); 

mat_files = arrayfun(@(x) fullfile(mat_path, x.name), dir(fullfile(mat_path, '*.mat')), 'uni', 0);

%% Configs 

stim_info = struct; % analyze all
stim_info.pulse_on_time         = 0; % ms 
stim_info.pulse_off_time        = Inf; % ms 

spike_params = struct; 
spike_params.peak_prom          = 10; 
spike_params.min_dist           = 0.5; 
spike_params.min_height         = -20; 
spike_params.min_width          = 0.05;  
spike_params.analyze_in_stim    = true; 

AP_params = struct; 
AP_params.max_tpre              = 3; 
AP_params.min_dV_dt             = 8;  
AP_params.scale_min_dV_dt       = false; 
% AP_params.min_dV_dt             = 0.01;  
% AP_params.scale_min_dV_dt       = true; 

AP_params.min_pts_cross_thres   = 2; 
AP_params.max_tpost             = 10; 

AP_params.rise_start_APfactor   = 0.1;
AP_params.rise_stop_APfactor    = 0.9;
AP_params.rise_interp_factor    = 100;
AP_params.fall_start_APfactor   = 0.9;
AP_params.fall_stop_APfactor    = 0.1; 
AP_params.fall_interp_factor    = 100;

configs.stim_info = stim_info;
configs.spike_params = spike_params;
configs.AP_params = AP_params; 

%% Extract
analysis_results = cellfun(@(x) process_data_per_file(x, configs, debug_fig_path), mat_files, 'uni', 0); 

save(fullfile(proc_path, 'analysis-summary.mat'), ...
    'analysis_results', '*path');


%% Helper function 
function summary_struct = process_data_per_file(file_name, config_struct, debug_fig_path)

if ~exist('debug_fig_path', 'var')
    debug_fig_path = nan;
end

load(file_name, 'recordings');
file_main_name = split(file_name, filesep);
file_main_name = file_main_name{end};

n_trace = length(recordings);
len_trace = length(recordings(1).data);
dt = recordings(1).dt * 1e3;
t = (0:len_trace-1)*dt;

init_src_struct = struct(...
    'configs', struct, ...
    'spike_properties', struct,...
    'AP_properies', struct); 

init_sum_struct = struct(...
    'Vthres_max', nan, ...
    'src', init_src_struct);

summary_struct = repmat(init_sum_struct, [n_trace,1]);


spike_params = config_struct.spike_params;
stim_info = config_struct.stim_info;
AP_params = config_struct.AP_params;

for i = 1:n_trace
    trace_struct = recordings(i);

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
    

    summary_struct(i).Vthres_max = max(AP_properies.AP_thres_Vm); 
    
    summary_struct(i).src.configs = config_struct;
    summary_struct(i).src.spike_properties = spike_properties;
    summary_struct(i).src.AP_properies = AP_properies;
    
    if ~isnan(debug_fig_path)
        
        trace_ttl = regexprep(file_main_name, '.mat', sprintf('-trace%02d', i));
        fig_name = fullfile(debug_fig_path, [trace_ttl '.pdf']);
        figure('units', 'normalized', 'position', [0,0,1,0.5], 'visible', 'off'); hold on;
        plot(t, Vm, '-k');
        plot(AP_properies.AP_time,AP_properies.AP_Vm, 'ro', 'MarkerFaceColor', 'r', 'displayname', 'spikes', 'tag', 'showlegend');
        plot(AP_properies.AP_thres_time,AP_properies.AP_thres_Vm, 'bo', 'MarkerFaceColor', 'b', 'displayname', 'thresholds',  'tag', 'showlegend');
        xlim([0, 0.6*t(end)]); ylim([-85,50]);
        xlabel('time (ms)'); ylabel('$V_m$ (mV)');
        title(regexprep(trace_ttl, '_', '-'));
        legend(findobj(gca, 'tag', 'showlegend'));
        despline;
        exportgraphics(gcf, fig_name);
        close;
    end
end
    

end
