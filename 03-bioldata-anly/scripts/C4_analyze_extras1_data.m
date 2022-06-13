clc; clear; close all;
run start_up.m

%% Define data paths
main_data_path = 'data/extra-s1/';
proc_path = fullfile(main_data_path, 'processed'); 

%%
load(fullfile(proc_path, 'analysis-summary.mat'), 'analysis_results');

%% 
rm_nan = @(x) x(~isnan(x));
get_all_Vthres = @(S) rm_nan(cell2mat(arrayfun(@(x) x.src.AP_properies.AP_thres_Vm, S, 'uni', 0)));
all_Vthres = cellfun(@(x) get_all_Vthres(x), analysis_results, 'uni', 0);

%% each cell has equal number of protocols, so just take mean 
Vthres_max = cellfun(@max, all_Vthres);

save(fullfile(proc_path, 'hansel-pham-extra-s1-vthres'), 'all_Vthres', 'Vthres_max')
