clc; clear;
run start_up.m;

%% Set constant options
N_Y = 3;
N_X = N_Y * 50;

num_Wrand = 10;

num_dthres = 50;
thres_baseline = 0.5;

num_inputs = 2000;
sigma_noise = 0.8;

rng_dthres = [-1,1];

% NAN as a holder of varied parameter options 
def_opts =  struct(...
    'N_X', N_X, ...
    'N_Y', N_Y, ...
    'num_Wrand', num_Wrand, ...
    'num_dthres', num_dthres, ...
    'thres_baseline', thres_baseline, ...
    'num_inputs', num_inputs, ...
    'sigma_noise', sigma_noise, ...
    'rng_dthres', rng_dthres, ...
    'alpha_W', nan, ...
    'percent_complete_input', nan, ...
    'num_overlap_per_group', nan);

%% Set paramater options

num_overlap_per_group_vec = 0:5:45;
percent_complete_input_vec = 0.1:0.1:1;
alpha_W_vec = linspace(0, 1, 50);

param_opts = struct; 
param_opts.num_overlap_per_group = num_overlap_per_group_vec;
param_opts.percent_complete_input = percent_complete_input_vec;
param_opts.alpha_W = alpha_W_vec;

%% Run
results_dim_description = 'dthres x alpha_W x percent_complete_input x num_overlap_per_group'; 
results_all = cell(1, length(alpha_W_vec), length(percent_complete_input_vec), length(num_overlap_per_group_vec));

delete(findall(0, 'Tag','progress_waitbar_run'));
wait_bar_handle = waitbar(0,'Please wait...','Tag','progress_waitbar_run');
set(get(findobj(wait_bar_handle,'Type','axes'),'Title'), ...
    'FontSize', 10, 'Interpreter', 'tex', 'FontName', 'FreeSans'); 
pause(0.5); 

%%

total_niters = length(num_overlap_per_group_vec) * length(percent_complete_input_vec);
current_niters = 1;

t0_run = tic;

for k = 1:length(num_overlap_per_group_vec)
    num_overlap_per_group = num_overlap_per_group_vec(k);
    
    for j = 1:length(percent_complete_input_vec)
        percent_complete_input = percent_complete_input_vec(j);
        
        parfor i = 1:length(alpha_W_vec)
            alpha_W = alpha_W_vec(i);
            
            opts = struct(...
                'N_X', N_X, ...
                'N_Y', N_Y, ...
                'num_Wrand', num_Wrand, ...
                'alpha_W', alpha_W, ...
                'num_dthres', num_dthres, ...
                'thres_baseline', thres_baseline, ...
                'num_inputs', num_inputs, ...
                'sigma_noise', sigma_noise, ...
                'rng_dthres', rng_dthres, ...
                'percent_complete_input', percent_complete_input, ...
                'num_overlap_per_group', num_overlap_per_group);
            
            results_all{1,i,j,k} = get_specperf(opts);
            
        end
        
        
        elapsed_run = toc(t0_run)/60;
        percent_run = current_niters/total_niters;
        est_remain = elapsed_run/percent_run - elapsed_run;
        
        waitbar(percent_run, wait_bar_handle, ...
            sprintf('Progress: %.0f%% (%d/%d) ..... Took %.1f min. ETA %.1f min', ...
            100*percent_run, current_niters, total_niters, ...
            elapsed_run, est_remain));
        
        current_niters = current_niters + 1;
        
    end
    
end

close(wait_bar_handle);

fprintf('Progress: 100%% (%d iters) ..... Took total %.1f min. \n', ...
    total_niters, elapsed_run);

results_all = cell2mat(results_all);

results_all = structarray_to_struct(results_all, 0);
results_all = structfun(@cell2mat, results_all, 'uni', 0);

save('data/recog_data.mat', ...
    'results_all', 'results_dim_description', ...
    'def_opts', 'param_opts'); 


