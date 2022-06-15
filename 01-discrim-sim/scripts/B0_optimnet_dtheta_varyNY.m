clc; clear; close all;
run startup.m

%% Define paths
data_path = 'data/vary10XNY';
fig_path = 'figures'; 
%% Testing 1 exampple

N_X = 10; 
N_Y = 6; 
dtheta_vec = linspace(-1, 1, 100);
n_iter = 3000; %5000;
eta_grad = 1e0; 
epsi_deltafun = 1e-2; 

results = run_1_optim_and_applydtheta(N_X, N_Y, dtheta_vec, n_iter, eta_grad, epsi_deltafun);

plot_one_example(results);

%% Define parameters for multiple simulations

% for testing par sim with smaller set
% N_X = 10; 
% N_Y_vec = [2,5,9,10]; 
% n_sim = 20;

N_X = 10; 
N_Y_vec = 2:20; 
n_sim = 500; 

dtheta_vec = linspace(-1, 1, 200);
n_iter = 5000;
eta_grad = 1e0; 
epsi_deltafun = 1e-2; 

num_N_Y = length(N_Y_vec);

%% Running multiple 

if isempty(gcp('nocreate'))
    parpool('local', 6);
end

results = cell(num_N_Y, 1); 

t0 = tic;
for i = 1:num_N_Y
    ti = tic; 
    N_Y = N_Y_vec(i); 
    fprintf('- Running simulations with N_Y=%d ... \t', N_Y);
    saved_filename = fullfile(data_path, sprintf('sim-NY=%d.mat', N_Y));
    results{i} = run_parmulti_optim_and_applydtheta(saved_filename, n_sim, ...
        N_X, N_Y, dtheta_vec, n_iter, eta_grad, epsi_deltafun);
    
    fprintf('took %.2f minutes.\n', toc(ti)/60);
end

fprintf('--> Total took %.2f minutes.\n', toc(t0)/60);

%% Plot an example from saved simulations
N_Y_sel = 5;
select_results = load(fullfile(data_path, sprintf('sim-NY=%d.mat', N_Y_sel)), 'simulations').simulations;
example_results = select_results{30};
plot_one_example(example_results);

exportgraphics(gcf, fullfile('figures', 'demo-optim.pdf'), 'ContentType', 'vector');

%% Load and plot entropies, num_unqs
data_files = arrayfun(@(x) fullfile(x.folder, x.name), dir(fullfile(data_path, '*.mat')), 'uni', 0); 
[~, sorted_by] = sort(cellfun(@(x) load(x).configurations.N_Y, data_files)); % sort by N_Y 
results = cellfun(@(x) load(x).analyses, data_files(sorted_by), 'uni', 0); 

results_agg = struct(...
    'init', {cellfun(@(x) x.init, results,'uni',0)}, ...
    'best', {cellfun(@(x) x.best, results,'uni',0)} ...
    );

prog_fields = fieldnames(results_agg);
prog_linestyles = struct('init', '-', 'best', '-');
        
z_sem = 2; 
fill_alpha = 0.2;
cmap = jet(num_N_Y)*0.9; 

figure; 
for i = 1:length(prog_fields)
    prog_field = prog_fields{i};
    pltsty = {'linestyle', prog_linestyles.(prog_field)};
    lgdn_opts = {true, '', {'NumColumns', 2, 'Location', 'northeast', 'FontSize', 12}};
    
    subplot(2,2,i); hold on;
    ttl = sprintf('Entropy of output patterns (%s)', prog_field);
    plot_single_result_panel(N_Y_vec, dtheta_vec, results_agg.(prog_field), 'entropy', ...
        z_sem, cmap, fill_alpha, pltsty, ttl, 'entropy (bits)', false);
    set(gca, 'tag', 'entropy');
    
    subplot(2,2,i+2); hold on;
    lgdn_opts_on = {{false}, lgdn_opts};
    ttl = sprintf('Number of unique output patterns (%s)', prog_field);
    plot_single_result_panel(N_Y_vec, dtheta_vec, results_agg.(prog_field), 'num_unq', ...
        z_sem, cmap, fill_alpha, pltsty, ttl, '\# patterns', lgdn_opts_on{(i==1) +1}{:})
    set(gca, 'tag', 'num_unq');
end

despline('all');

linkaxes(findall(gcf,'type','axes'),'x');
cellfun(@(x) linkaxes(findall(gcf, 'type', 'axes', 'tag', x), 'y'), {'entropy', 'num_unq'})

pause(1);
exportgraphics(gcf, fullfile(fig_path, 'discrim-ent-unq-10X-varyNY.pdf'))

%% Load data to analyze further
num_N_Y = length(data_files);

N_Y_vec = zeros(num_N_Y,1);
results_agg = cell(num_N_Y, 1); 

sm_wins = struct('entropy', 5, 'num_unq', 10);
net_fields = {'theta', 'W',}';

tic
for i = 1:num_N_Y
    dat = load(data_files{i}); 
    res = struct();
    
    N_Y_vec(i) = dat.configurations.N_Y;
    dtheta_vec = dat.configurations.dtheta_vec;
    
    L1_fields = fieldnames(dat.simulations{1}.net_config);
    L2_fields = fieldnames(dat.simulations{1}.apply_dtheta.(L1_fields{1}));
    
    for j1 = 1:length(L1_fields)
        F1 = L1_fields{j1}; 
        res.(F1) = struct;
        
        for j2 = 1:length(net_fields)
            nF = net_fields{j2};
            res.(F1).net.(nF) = cellfun(@(x) mean(x.net_config.(F1).(nF), 'all'), dat.simulations);
        end
        
        for j2 = 1:length(L2_fields)
            F2 = L2_fields{j2};
            res.(F1).(F2) = structarray_to_struct(...
                cellfun(@(D) analyze_dthetacurve(dtheta_vec, ...
                    D.apply_dtheta.(F1).(F2), sm_wins.(F2)), ...
                    dat.simulations,  'uni', 1));
        end
    end
    
    results_agg{i} = res;
end

toc

[N_Y_vec, sorted_inds_NY] = sort(N_Y_vec);
results_agg = results_agg(sorted_inds_NY);
results_agg = structarray_to_struct(vertcat(results_agg{:}));
results_agg = structfun(@(S) structarray_to_struct(S, 0), results_agg, 'uni', 0);

prog_fields = L1_fields;
anly_fields = [{'net'}; L2_fields];

for i1 = 1:length(prog_fields)
    F1 = prog_fields{i1};

    for i2 = 1:length(anly_fields)
        F2 = anly_fields{i2};
        results_agg.(F1).(F2) = structfun(@(x) horzcat(x{:})', ...
            structarray_to_struct([results_agg.(F1).(F2){:}], 0), 'uni',0);

    end
end

%% Plot
range_NY_plt = 2:20;
ind_NY_plt = arrayfun(@(x) find_nearest(N_Y_vec, x, 'ind'), range_NY_plt);
cmap = jet(length(range_NY_plt))*0.9; 

plt_pairs = {...
    struct('pairs', {{'best-net-theta', 'best-net-W'}}); ...
    struct('pairs', {{'init-net-theta', 'best-net-theta'}}, 'idline', true); ...
    struct('pairs', {{'init-net-W', 'best-net-W'}}, 'idline', true); ...
    struct('pairs', {{'best-entropy-dtheta_at_max', 'best-entropy-max'}}); ...
    struct('pairs', {{'best-net-theta', 'best-entropy-dtheta_at_max'}}); ...
    struct('pairs', {{'best-entropy-mean_at_negdtheta', 'best-entropy-mean_at_posdtheta'}}, 'idline', true); ...
    struct('pairs', {{'best-num_unq-dtheta_at_max', 'best-num_unq-max'}}); ...
    struct('pairs', {{'best-net-theta', 'best-num_unq-dtheta_at_max'}}); ...
    struct('pairs', {{'best-num_unq-mean_at_negdtheta', 'best-num_unq-mean_at_posdtheta'}}, 'idline', true); ...
    };

fields_to_labels = struct(...
    'best', 'best', ...
    'init', 'init', ...
    'net', '', ...
    'W', '$W$', ...
    'theta', '$\theta$', ...
    'entropy', 'ent', ...
    'num_unq', 'unq', ...
    'max', 'max', ...
    'dtheta_at_max', '$\Delta\theta_{\mathrm{max}}$', ...
    'mean', 'mean', ...
    'mean_at_negdtheta', 'mean$_{\Delta\theta < 0}$', ...
    'mean_at_posdtheta', 'mean$_{\Delta\theta > 0}$' ...
    );

parse_labels_in_fields = @(x) strsplit(x,'-');
parse_fields_to_labels = @(x) strjoin(cellfun(@(k) fields_to_labels.(k), x, 'uni', 0), '-');

num_pairs = length(plt_pairs);
nrows = ceil(sqrt(num_pairs)); ncols = ceil(num_pairs/nrows); 

figure; 
for i = 1:num_pairs
    subplot(nrows, ncols, i); hold on;
    pair_ith = plt_pairs{i}.pairs;
    
    idline = false; 
    if isfield(plt_pairs{i}, 'idline'), idline = plt_pairs{i}.idline; end
    
    xlbl = pair_ith{1}; xfields = parse_labels_in_fields(xlbl);
    ylbl = pair_ith{2}; yfields = parse_labels_in_fields(ylbl);
    X_mat = getfield(results_agg, xfields{:});
    Y_mat = getfield(results_agg, yfields{:});
    X_mat = X_mat(ind_NY_plt,:);
    Y_mat = Y_mat(ind_NY_plt,:);
    if idline 
        rng_all_mats = [min(min(X_mat(:)), min(Y_mat(:))), max(max(X_mat(:)), max(Y_mat(:)))];
        plot(rng_all_mats, rng_all_mats, '-', 'linewidth', 1.5, 'color', [0.1,0.1,0.1,0.8]);
        
        daspect([1,1,1]);
    end
    
    scatter_with_colorgrads(X_mat, Y_mat, cmap, 40, 0.02)
    xlabel(parse_fields_to_labels(xfields)); 
    ylabel(parse_fields_to_labels(yfields));

end

colormap(cmap);
cbar = colorbar('eastoutside');
cbar.Position = cbar.Position .* [1,1,1,3] + [0.1,0.1,0,0];
title(cbar, '$N_Y$', 'interpreter', 'latex')
caxis([min(range_NY_plt), max(range_NY_plt)]);
set(cbar, 'ticks', range_NY_plt);

despline('all');
pause(1);

exportgraphics(gcf, fullfile(fig_path, 'discrim-pairs-anly-10X-varyNY.pdf'))

%% Helper plot functions
function scatter_with_colorgrads(X_mat, Y_mat, cmap, markersize, markeralpha)
    if ~exist('markersize', 'var'), markersize = 50; end 
    if ~exist('markeralpha', 'var'), markeralpha = 0.1; end
    mean_mkalp = max(markeralpha*5, 0.4);
    mean_mksz = round(markersize*1.5);
    
    arrayfun(@(i) scatter(X_mat(i,:), Y_mat(i,:), ...
        markersize, cmap(i,:), 'filled', 'o', ...
        'MarkerEdgeAlpha', 0, 'MarkerFaceAlpha', markeralpha),  1:size(X_mat,1));
    
    arrayfun(@(i) scatter(mean(X_mat(i,:)), mean(Y_mat(i,:)), ...
        mean_mksz, cmap(i,:), 'filled', 'o', ...
        'markeredgecolor', 'k', 'LineWidth', 0.8, ...
        'MarkerEdgeAlpha', 0.5, 'MarkerFaceAlpha', mean_mkalp),  1:size(X_mat,1));

end


