clc; clear; close all; 
run startup.m

%% Define path 
% variation of N_Y and k
data_path = 'data/varyNYk'; 

if ~exist(data_path, 'dir')
    mkdir(data_path);
end

%% Define parameters
N_X = 10;
N_Y_vec = [6, 10, 14]; 
k_vec = 2:8;

n_dtheta = 100;
n_sim = 100; 
n_iter = 2500; 

%% Simulation initialization
num_NY = length(N_Y_vec);
num_k = length(k_vec);
dtheta_vec = linspace(-1, 1, n_dtheta); 

Xids = 0:(2^N_X-1);
Xs = de2bi(Xids)';

%% Run simulation 
analyses = cell(num_NY, num_k, n_dtheta, n_sim);
analyses(:) = {struct()}; 

for i_ny = 1:num_NY
    N_Y = N_Y_vec(i_ny);
    tic;
    
    ppm = ParforProgMon(sprintf('N_Y=%d', N_Y), n_sim);
    % run multiple sims
    parfor i_sim = 1:n_sim
        
        % optim
        net_config = run_1_optim(N_X, N_Y, n_iter);
        net_fields = {'init', 'best'};
        
        % either init or best 
        for i_net = 1:length(net_fields)
            nf = net_fields{i_net};
            
            % get net config
            W = net_config.(nf).W;
            theta = net_config.(nf).theta;
            
            % apply dtheta
            for i_dtheta = 1:n_dtheta
                d_theta = dtheta_vec(i_dtheta);
                theta_app = theta + d_theta;
                
                Ys = W * Xs > theta_app;
                
                for i_k = 1:num_k
                    k = k_vec(i_k);
                    sub_X = sum(Xs(1:k,:), 1) == k;
                    Y_subX = Ys(:,sub_X)';
                    
                    Jd = get_binary_pairwise_measure_distribution(Y_subX, @calculate_Jaccard_distance);
                    Hd = get_binary_pairwise_measure_distribution(Y_subX, @calculate_Hamming_distance)/N_Y;
                    
                    tmp = struct(...
                        'k', k, ...
                        'Y_mean', mean(Y_subX,'all'), ...
                        'Jd_mean', mean(Jd), ...
                        'Jd_med', median(Jd), ...
                        'Jd_sd', std(Jd), ...
                        'Hd_mean', mean(Hd), ...
                        'Hd_med', median(Hd), ...
                        'Hd_sd', std(Hd),...
                        'Jd_mean_nonan', mean(Jd, 'omitnan'), ...
                        'Jd_med_nonan', median(Jd, 'omitnan'), ...
                        'Jd_sd_nonan', std(Jd, 'omitnan') ...
                        );
                    
                    analyses{i_ny, i_k, i_dtheta, i_sim}.(nf) = tmp; 
                end
                
                
            end
        end
        
        ppm.increment();
    end
    toc;
    delete(ppm);
end

%% Save results 
config = struct(...
    'N_X', N_X, ...
    'N_Y_vec', N_Y_vec, ...
    'k_vec', k_vec, ...
    'n_dtheta', n_dtheta, ...
    'n_sim', n_sim, ...
    'n_iter', n_iter, ...
    'dtheta_vec', dtheta_vec);

save(fullfile(data_path, 'sim_results.mat'), 'analyses', 'config');

%% Summary analyses (take means of trials)
load(data_file, 'analyses', 'config');
results = cell2mat(analyses); clear analyses; 
results = structarray_to_struct(results);
results = structfun(@structarray_to_struct, results, 'uni', 0);
summary_analyses = structfun(@(s) structfun(@(x) mean(x,4), s, 'uni', 0), results, 'uni', 0);

save(fullfile(data_path, 'summary_analyses.mat'), 'config', 'summary_analyses'); 