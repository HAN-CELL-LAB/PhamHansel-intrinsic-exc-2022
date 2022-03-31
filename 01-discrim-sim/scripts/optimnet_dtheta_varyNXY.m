clc; clear; close all; 

addpath(genpath('extpkg'));
addpath(genpath('functions'));
graphic_setdefault(15, ...
    'DefaultAxesMinorGridAlpha', 0.05, ...
    'DefaultAxesMinorGridLineStyle', '-', ...
    'DefaultTextInterpreter', 'latex', ...
    'DefaultLegendInterpreter', 'latex', ...
    'DefaultStemMarkerSize', 1, ...
    'DefaultStemlineWidth', 1.5, ...
    'DefaultFigureWindowStyle','normal');

%% Paths
data_path = 'data/varyNXY'; 
fig_path = 'figures/varyNXY';

if ~exist(data_path, 'dir')
    mkdir(data_path);
end

if ~exist(fig_path, 'dir')
    mkdir(fig_path);
end

%% Testing 1 exampple
% N_X = 15; 
% N_Y = 10; 
% dtheta_vec = linspace(-1, 1, 100);
% n_iter = 3000; %5000;
% eta_grad = 1e0; 
% epsi_deltafun = 1e-2; 
% 
% tic
% results = run_1_optim_and_applydtheta(N_X, N_Y, dtheta_vec, n_iter, eta_grad, epsi_deltafun);
% toc
% plot_one_example(results);

%% Combinations
% params_to_vary = {...
%     struct('label', 'N_Y', 'vec', [3, 5, 8, 10, 12, 15, 17, 20]), ...
%     struct('label', 'N_X', 'vec', [5, 10, 15]) ...
% };
% 
% [combinations, labels, comb_indices] = return_combination(params_to_vary{:}); 
% comb_table = array2table(combinations, 'VariableNames', labels);
% 
% n_comb = size(comb_table, 1);
% 
% n_sim = 100; 
% 
% dtheta_vec = linspace(-1, 1, 200);
% 
% n_iter = 2500;
% eta_grad = 1e0; 
% epsi_deltafun = 1e-2;  
% 
% params_config = struct(...
%     'params_to_vary', {params_to_vary}, ...
%     'combination_table', comb_table, ...
%     'n_comb', n_comb, ...
%     'n_sim', n_sim, ...
%     'dtheta_vec', dtheta_vec, ...
%     'n_iter', n_iter, ...
%     'eta_grad', eta_grad, ...
%     'epsi_deltafun', epsi_deltafun);
% 
% save(fullfile(data_path, 'params_config.mat'), 'params_config');

%% Running multiple 

% if isempty(gcp('nocreate'))
%     parpool('local', 6);
% end
% 
% 
% t0 = tic;
% for i_comb = 1:n_comb
%     ti = tic;
%     
%     comb_ith = comb_table(i_comb,:);
%     N_X = comb_ith.N_X;
%     N_Y = comb_ith.N_Y;
%     fprintf('- (%02d) Running simulations with: N_X=%d \t N_Y=%d ... \t', ...
%         i_comb, N_X, N_Y);
%     
%     saved_filename = fullfile(data_path, sprintf('sim_%03d.mat', i_comb));
%     run_parmulti_optim_and_applydtheta(saved_filename, n_sim, ...
%         N_X, N_Y, dtheta_vec, n_iter, eta_grad, epsi_deltafun);
%     
%     fprintf('took %.2f minutes.\n', toc(ti)/60);
% end
% 
% fprintf('--> Total took %.2f minutes.\n', toc(t0)/60);

%% Load to plot

load(fullfile(data_path, 'params_config.mat'), 'params_config');
combination_table = params_config.combination_table;
dtheta_vec = params_config.dtheta_vec;
num_comb = size(combination_table,1);
N_X_vec = unique(combination_table.N_X);
N_Y_vec = unique(combination_table.N_Y);
num_N_X = length(N_X_vec);
num_N_Y = length(N_Y_vec);

%%
all_analyses = cell(num_comb, 1); 
for i_comb = 1:num_comb
    comb_ith = combination_table(i_comb, :);
    dat_file = fullfile(data_path, sprintf('sim_%03d.mat', i_comb));
    data_ith = load(dat_file, 'configurations', 'analyses');
    
    if data_ith.configurations.N_X ~= comb_ith.N_X || ...
            data_ith.configurations.N_Y ~= comb_ith.N_Y
        error('Mismatch config'); 
    end
    
    all_analyses{i_comb} = data_ith.analyses;
end

analyses_aggr = struct(...
    'init', {cellfun(@(x) x.init, all_analyses,'uni',0)}, ...
    'best', {cellfun(@(x) x.best, all_analyses,'uni',0)} ...
    );

%%
prog_fields = fieldnames(analyses_aggr);
num_progs = length(prog_fields);
prog_linestyles = struct('init', '-', 'best', '-');

z_sem = 2; 
fill_alpha = 0.2;
cmap = jet(num_N_Y)*0.9; 

nrows = num_progs; 
ncols = num_N_X; 

figure; 
annotation('textbox', 'string', 'Entropy of output patterns', ...
    'units', 'normalized', 'position', [0,0.97,1,0.03], ...
    'LineStyle', 'none', 'FontSize', 20, 'Interpreter', 'latex', ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'cap');

for i_nx = 1:num_N_X
    N_X_sel = N_X_vec(i_nx); 
    sel_ind_N_X = find(combination_table.N_X == N_X_sel); 
    N_Y_sel = combination_table.N_Y(sel_ind_N_X); 
    
    [N_Y_sel, sorted_ind] = sort(N_Y_sel);
    
    for i_p = 1:num_progs
        prog_field = prog_fields{i_p};
        
        anly_sel = analyses_aggr.(prog_field)(sel_ind_N_X);
        anly_sel = anly_sel(sorted_ind);
        
        pltsty = {'linestyle', prog_linestyles.(prog_field)};
        lgdn_opts = {true, '', {'NumColumns', 1, 'Location', 'northeast', 'FontSize', 12}};
        
        splt_cnt = i_nx + (i_p-1)*ncols;
        subplot(nrows, ncols, splt_cnt); hold on;
        
        lgdn_opts_on = {{false}, lgdn_opts};
        
        anly_field = 'entropy';
        ylbl = 'entropy (bits)';
        
        if i_nx == 1
            ylbl = {sprintf('\\textit{%s}', prog_field), ylbl}; %#ok<AGROW>
        end
        
        ttl = sprintf('$N_X=%d$', N_X_sel);
        
        plot_single_result_panel(N_Y_sel, dtheta_vec, anly_sel, anly_field, ...
            z_sem, cmap, fill_alpha, pltsty, ttl, ylbl,  ...
            lgdn_opts_on{(splt_cnt==nrows*ncols) +1}{:})
        set(gca, 'tag', anly_field);
        
    end
end

despline('all');

linkaxes(findall(gcf,'type','axes'),'xy');

pause(1);
exportgraphics(gcf, fullfile(fig_path, 'anly-results-entropy.pdf'), 'ContentType', 'vector')

%%
figure; 
annotation('textbox', 'string', 'Entropy of output patterns', ...
    'units', 'normalized', 'position', [0,0.97,1,0.03], ...
    'LineStyle', 'none', 'FontSize', 20, 'Interpreter', 'latex', ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'cap');

for i_nx = 1:num_N_X
    N_X_sel = N_X_vec(i_nx); 
    sel_ind_N_X = find(combination_table.N_X == N_X_sel); 
    N_Y_sel = combination_table.N_Y(sel_ind_N_X); 
    
    [N_Y_sel, sorted_ind] = sort(N_Y_sel);
    
    for i_p = 1:num_progs
        prog_field = prog_fields{i_p};
        
        anly_sel = analyses_aggr.(prog_field)(sel_ind_N_X);
        anly_sel = anly_sel(sorted_ind);
        
        pltsty = {'linestyle', prog_linestyles.(prog_field)};
        lgdn_opts = {true, '', {'NumColumns', 1, 'Location', 'northeast', 'FontSize', 12}};
        
        splt_cnt = i_nx + (i_p-1)*ncols;
        subplot(nrows, ncols, splt_cnt); hold on;
        
        lgdn_opts_on = {{false}, lgdn_opts};
        
        anly_field = 'num_unq';
        ylbl = '\# patterns';
        
        if i_nx == 1
            ylbl = {sprintf('\\textit{%s}', prog_field), ylbl}; %#ok<AGROW>
        end
        
        ttl = sprintf('$N_X=%d$', N_X_sel);
        
        plot_single_result_panel(N_Y_sel, dtheta_vec, anly_sel, anly_field, ...
            z_sem, cmap, fill_alpha, pltsty, ttl, ylbl,  ...
            lgdn_opts_on{(splt_cnt==nrows*ncols) +1}{:})
        set(gca, 'tag', anly_field);
        
    end
end

despline('all');

linkaxes(findall(gcf,'type','axes'),'xy');

pause(1);
exportgraphics(gcf, fullfile(fig_path, 'anly-results-num_unq.pdf'), 'ContentType', 'vector')


%% Plot an example from saved simulations
% data_path = 'data'; 
% N_Y_sel = 5;
% select_results = load(fullfile(data_path, sprintf('sim-NY=%d.mat', N_Y_sel)), 'simulations').simulations;
% num_sim = length(select_results); 
% example_results = select_results{randi(num_sim)};
% plot_one_example(example_results);
% 
% % exportgraphics(gcf, fullfile('figures', 'example-sim01.pdf'), 'ContentType', 'vector')

%% Load data to analyze further
for i_nx = 1:num_N_X
    N_X_sel = N_X_vec(i_nx); 
    sel_ind_N_X = find(combination_table.N_X == N_X_sel); 
    N_Y_sel = combination_table.N_Y(sel_ind_N_X); 
    
    data_files = arrayfun(@(x) fullfile(data_path, sprintf('sim_%03d.mat', x)), sel_ind_N_X, 'uni', 0);
    
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
range_NY_plt = N_Y_vec;
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
    
    scatter_with_colorgrads(X_mat, Y_mat, cmap, 50, 0.1)
    if i == 1
        title(sprintf('$N_X=%d$', N_X_sel));
    end
    
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
% 
% % savefig(gcf, fullfile('figures', 'pairs-results.fig')); 
% 
exportgraphics(gcf, fullfile(fig_path, sprintf('pairs-results-all_NX=%d.pdf', N_X_sel)), 'ContentType', 'image')

end
%%
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


