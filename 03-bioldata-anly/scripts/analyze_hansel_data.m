clc; clear; close all;
run start_up.m

%% Define data paths
main_data_path = 'data/hansel-data/';
proc_path = fullfile(main_data_path, 'processed'); 

%%
load(fullfile(proc_path, 'analysis-summary.mat'), 'analysis_results');

%% 
rm_nan = @(x) x(~isnan(x));
get_all_Vthres = @(S) rm_nan(cell2mat(arrayfun(@(x) x.src.AP_properies.AP_thres_Vm, S, 'uni', 0)));
all_Vthres = cellfun(@(x) get_all_Vthres(x), analysis_results, 'uni', 0);

%%
figure; hold on; 
histogram(cell2mat(all_Vthres), 100)