clc; clear; close all; 
run start_up.m

%% Define data paths
main_data_path = 'data/gill-data/Layer II.III S1 BC IP'; 
INFO_path = fullfile(main_data_path, 'info');

load(fullfile(INFO_path, 'info-summary.mat'));
INFO_file = fullfile(INFO_path, 'info-summary.xlsx');

%% There are some missing so wont be analyze from this 
passive_table = cellfun(@(x) readtable(INFO_file, 'Sheet', x, 'PreserveVariableNames', true), {'Vm', 'Rin'}, 'uni', 0);