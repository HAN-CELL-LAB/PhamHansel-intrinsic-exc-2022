clc; clear; close all; 

addpath(genpath('extpkg'));
addpath(genpath('functions'));

num_workers = 40;
parcluster_obj = parcluster;

if isempty(gcp('nocreate'))
    parpool(parcluster_obj, num_workers);
end

section_sep = repmat('-', 1, 50);
fprintf('Information about cluster: \n');
disp(parcluster_obj); 
fprintf('%s\n', section_sep);

%% Paths
data_path = 'data/varyNXY-sim01'; 
fig_path = 'figures/varyNXY-sim01';

if ~exist(fig_path, 'dir')
    mkdir(fig_path);
end

%% Combinations
params_to_vary = {...
    struct('label', 'N_Y', 'vec', 2:2:20), ...
    struct('label', 'N_X', 'vec', [5, 6, 7, 8, 9, 10, 12, 14, 16]) ...
};

[combinations, labels, comb_indices] = return_combination(params_to_vary{:}); 
comb_table = array2table(combinations, 'VariableNames', labels);

n_comb = size(comb_table, 1);

n_sim = 400; 

dtheta_vec = linspace(-1, 1, 200);

n_iter = 2500;
eta_grad = 1e0; 
epsi_deltafun = 1e-2;  

params_config = struct(...
    'params_to_vary', {params_to_vary}, ...
    'combination_table', comb_table, ...
    'n_comb', n_comb, ...
    'n_sim', n_sim, ...
    'dtheta_vec', dtheta_vec, ...
    'n_iter', n_iter, ...
    'eta_grad', eta_grad, ...
    'epsi_deltafun', epsi_deltafun);

save(fullfile(data_path, 'params_config.mat'), 'params_config');

%% Running multiple 
    
sim_progress = comb_table;
sim_progress.elapsedminutes = zeros(n_comb,1);
sim_progress.datafile = repmat("",n_comb,1);

t0 = tic;

fprintf('\n=============== START SIMULATION ===============\n');
for i_comb = 1:n_comb
    ti = tic;
    
    comb_ith = comb_table(i_comb,:);
    N_X = comb_ith.N_X;
    N_Y = comb_ith.N_Y;
    fprintf('- (%02d) Running simulations with: N_X=%02d    N_Y=%02d ... \t', ...
        i_comb, N_X, N_Y);
    
    saved_filename = fullfile(data_path, sprintf('sim_%03d.mat', i_comb));
    
    run_parmulti_optim_and_applydtheta(saved_filename, n_sim, ...
        N_X, N_Y, dtheta_vec, n_iter, eta_grad, epsi_deltafun);
    
    fprintf('took %.2f minutes.\n', toc(ti)/60);
    
    sim_progress.elapsedminutes(i_comb) = toc(ti)/60;
    sim_progress.datafile(i_comb) = string(saved_filename);
    
    save(fullfile(data_path, 'sim_progress.mat'), 'sim_progress');
end

fprintf('--> Total took %.2f minutes.\n', toc(t0)/60);


fprintf('\n=============== FINISH SIMULATION ===============\n');