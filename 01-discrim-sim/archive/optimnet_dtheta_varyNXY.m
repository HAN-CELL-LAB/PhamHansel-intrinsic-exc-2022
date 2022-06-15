clc; clear; close all; 
run startup.m

%% Paths
data_path = 'data/varyNXY'; 
fig_path = 'figures';

if ~exist(data_path, 'dir')
    mkdir(data_path);
end

%% Testing 1 exampple
N_X = 15; 
N_Y = 10; 
dtheta_vec = linspace(-1, 1, 100);
n_iter = 3000; %5000;
eta_grad = 1e0; 
epsi_deltafun = 1e-2; 

tic
results = run_1_optim_and_applydtheta(N_X, N_Y, dtheta_vec, n_iter, eta_grad, epsi_deltafun);
toc
plot_one_example(results);

%% Combinations
params_to_vary = {...
    struct('label', 'N_Y', 'vec', [3, 5, 8, 10, 12, 15, 17, 20]), ...
    struct('label', 'N_X', 'vec', [5, 10, 15]) ...
};

[combinations, labels, comb_indices] = return_combination(params_to_vary{:}); 
comb_table = array2table(combinations, 'VariableNames', labels);

n_comb = size(comb_table, 1);

n_sim = 100; 

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

if isempty(gcp('nocreate'))
    parpool('local', 6);
end


t0 = tic;
for i_comb = 1:n_comb
    ti = tic;
    
    comb_ith = comb_table(i_comb,:);
    N_X = comb_ith.N_X;
    N_Y = comb_ith.N_Y;
    fprintf('- (%02d) Running simulations with: N_X=%d \t N_Y=%d ... \t', ...
        i_comb, N_X, N_Y);
    
    saved_filename = fullfile(data_path, sprintf('sim_%03d.mat', i_comb));
    run_parmulti_optim_and_applydtheta(saved_filename, n_sim, ...
        N_X, N_Y, dtheta_vec, n_iter, eta_grad, epsi_deltafun);
    
    fprintf('took %.2f minutes.\n', toc(ti)/60);
end

fprintf('--> Total took %.2f minutes.\n', toc(t0)/60);

%% Load to plot

load(fullfile(data_path, 'params_config.mat'), 'params_config');
combination_table = params_config.combination_table;
dtheta_vec = params_config.dtheta_vec;
num_comb = size(combination_table,1);
N_X_vec = unique(combination_table.N_X);
N_Y_vec = unique(combination_table.N_Y);
num_N_X = length(N_X_vec);
num_N_Y = length(N_Y_vec);

%% Load analyses 

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

%% Plot entropy
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
exportgraphics(gcf, fullfile(fig_path, 'discrim-ent-varyNXY.pdf'))

%% Plot number of unique patterns

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
        set(gca, 'tag', anly_field, 'yscale', 'log');
        
    end
end

despline('all');

linkaxes(findall(gcf,'type','axes'),'xy');

pause(1);
exportgraphics(gcf, fullfile(fig_path, 'discrim-unq-varyNXY.pdf'))


