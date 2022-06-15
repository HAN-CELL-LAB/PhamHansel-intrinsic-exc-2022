clc; clear; close all; 
run startup.m

%% Paths
data_path = 'data/varyNXY'; 
fig_path = 'figures/';

if ~exist(fig_path, 'dir')
    mkdir(fig_path);
end

%% Load to plot
load(fullfile(data_path, 'params_config.mat'), 'params_config');
combination_table = params_config.combination_table;
dtheta_vec = params_config.dtheta_vec;
num_comb = size(combination_table,1);
N_X_vec = unique(combination_table.N_X);
N_Y_vec = unique(combination_table.N_Y);
num_N_X = length(N_X_vec);
num_N_Y = length(N_Y_vec);

%% Aggregate analysis

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

%% Plot entropies of pattern distributions

prog_fields = fieldnames(analyses_aggr);
num_progs = length(prog_fields);
prog_linestyles = struct('init', '-', 'best', '-');

z_sem = 3; 
fill_alpha = 0.4;
cmap = jet(num_N_Y)*0.9; 


figure; 

N_X_sel_vec = [6, 10, 14]; 
N_X_sel_ind = arrayfun(@(nx) find_nearest(N_X_vec, nx, 'ind'), N_X_sel_vec);

nrows = num_progs; 
ncols = length(N_X_sel_ind);

for i_nx = 1:length(N_X_sel_ind)
    N_X_sel = N_X_sel_vec(i_nx); 
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
        
        xlabel('$\Delta \theta$');
        title(sprintf('$N_X=%d$', N_X_sel));
        set(gca, 'tag', anly_field);
        
    end
end

despline('all');

linkaxes(findall(gcf,'type','axes'),'xy');

exportgraphics(gcf, fullfile(fig_path, 'discrim-ent-varyNXY.pdf'))

%% Plot number of unique patterns 
% why, yes I am very lazy 

prog_fields = fieldnames(analyses_aggr);
num_progs = length(prog_fields);
prog_linestyles = struct('init', '-', 'best', '-');

z_sem = 3; 
fill_alpha = 0.4;
cmap = jet(num_N_Y)*0.9; 

figure; 

N_X_sel_vec = [6, 10, 14]; 
N_X_sel_ind = arrayfun(@(nx) find_nearest(N_X_vec, nx, 'ind'), N_X_sel_vec);

nrows = num_progs; 
ncols = length(N_X_sel_ind);
for i_nx = 1:length(N_X_sel_ind)
    N_X_sel = N_X_sel_vec(i_nx); 
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
        ylbl = '\# unq pats';
        
        if i_nx == 1
            ylbl = {sprintf('\\textit{%s}', prog_field), ylbl}; %#ok<AGROW>
        end
        
        ttl = sprintf('$N_X=%d$', N_X_sel);
        
        plot_single_result_panel(N_Y_sel, dtheta_vec, anly_sel, anly_field, ...
            z_sem, cmap, fill_alpha, pltsty, ttl, ylbl,  ...
            lgdn_opts_on{(splt_cnt==nrows*ncols) +1}{:})
        
        xlabel('$\Delta \theta$');
        title(sprintf('$N_X=%d$', N_X_sel));        
        set(gca, 'tag', anly_field, 'yscale', 'log');
        
    end
end

despline('all');

linkaxes(findall(gcf,'type','axes'),'xy');

exportgraphics(gcf, fullfile(fig_path, 'discrim-unq-varyNXY.pdf'))
