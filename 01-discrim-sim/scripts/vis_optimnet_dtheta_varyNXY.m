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
data_path = 'data/varyNXY-sim01'; 
fig_path = 'figures/varyNXY-sim01';

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

%% 

prog_fields = fieldnames(analyses_aggr);
num_progs = length(prog_fields);
prog_linestyles = struct('init', '-', 'best', '-');

z_sem = 3; 
fill_alpha = 0.4;
cmap = jet(num_N_Y)*0.9; 

nrows = num_progs; 
ncols = num_N_X; 

figure; 
annotation('textbox', 'string', 'Norm (ind) Entropy of output patterns', ...
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
%         ylbl = '';
        if i_nx == 1
%             ylbl = 'entropy (bits)';
            ylbl = {sprintf('\\textit{%s}', prog_field), ylbl}; %#ok<AGROW>
        end
        
        ttl = sprintf('$N_X=%d$', N_X_sel);
        
%         plot_single_result_panel(N_Y_sel, dtheta_vec, anly_sel, anly_field, ...
%             z_sem, cmap, fill_alpha, pltsty, ttl, ylbl,  ...
%             lgdn_opts_on{(splt_cnt==nrows*ncols) +1}{:})
        
        arrayfun(@(i_ny) plot(dtheta_vec, anly_sel{i_ny}.entropy.mean ./  log2(anly_sel{i_ny}.num_unq.mean), ...
            'color', cmap(i_ny,:)), 1:length(N_Y_sel));
        
%         arrayfun(@(i_ny) plot(dtheta_vec, anly_sel{i_ny}.entropy.mean ./  N_X_sel, ...
%             'color', cmap(i_ny,:)), 1:length(N_Y_sel));
%         
        xlabel('$\Delta \theta$');
        if i_nx == 1
%         ylabel({sprintf('\\textit{%s}', prog_field), 'ent/$N_X$'});
        
        ylabel({sprintf('\\textit{%s}', prog_field), 'ent/$\log_2(K)$'});
        end
        title(sprintf('$N_X=%d$', N_X_sel));
        
        set(gca, 'tag', anly_field);
        
    end
end

despline('all');

linkaxes(findall(gcf,'type','axes'),'xy');
% waitforbuttonpress;
% pause(1);
% exportgraphics(gcf, fullfile(fig_path, 'anly-results-entropy.pdf'), 'ContentType', 'vector')
% exportgraphics(gcf, fullfile(fig_path, 'anly-results-entropy.pdf'), 'ContentType', 'vector')

%%
prog_fields = fieldnames(analyses_aggr);
num_progs = length(prog_fields);
prog_linestyles = struct('init', '-', 'best', '-');

z_sem = 3; 
fill_alpha = 0.4;
cmap = jet(num_N_Y)*0.9; 

nrows = num_progs; 
ncols = num_N_X; 

figure; 
annotation('textbox', 'string', 'Norm Number of unique output patterns', ...
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
        ylbl = '\# unq pats';
        
        ylbl = '';
        if i_nx == 1
            
        ylbl = '\# unq pats';
            ylbl = {sprintf('\\textit{%s}', prog_field), ylbl}; %#ok<AGROW>
        end
        
        ttl = sprintf('$N_X=%d$', N_X_sel);
        
%         plot_single_result_panel(N_Y_sel, dtheta_vec, anly_sel, anly_field, ...
%             z_sem, cmap, fill_alpha, pltsty, ttl, ylbl,  ...
%             lgdn_opts_on{(splt_cnt==nrows*ncols) +1}{:})
        
        
        arrayfun(@(i_ny) plot(dtheta_vec, anly_sel{i_ny}.num_unq.mean ./  2^N_X_sel, ...
            'color', cmap(i_ny,:)), 1:length(N_Y_sel));
        
        xlabel('$\Delta \theta$');
        if i_nx == 1
        ylabel({sprintf('\\textit{%s}', prog_field), 'K/$2^{N_X}$'});
        end
        title(sprintf('$N_X=%d$', N_X_sel));
        


        
        set(gca, 'tag', anly_field);
        
    end
end

despline('all');

linkaxes(findall(gcf,'type','axes'),'xy');

%%
cmap = jet(num_N_Y)*0.9; 

nrows = num_progs;
ncols = 1; 

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
        
        splt_cnt = i_p; 
        subplot(nrows, ncols, splt_cnt); hold on;
        
        lgdn_opts_on = {{false}, lgdn_opts};
        
        anly_field = 'entropy';
        ylbl = 'entropy (bits)';
        
        ylbl = {sprintf('\\textit{%s}', prog_field), ylbl}; %#ok<AGROW>        
        ttl = '';
        
        cmap_adj = bound_minmax(cmap + 0.3*(1-i_nx/num_N_X), 0, 0.96);
        arrayfun(@(i_ny) plot(dtheta_vec, anly_sel{i_ny}.entropy.mean ./  log2(anly_sel{i_ny}.num_unq.mean), ...
            'color', cmap_adj(i_ny,:)), 1:length(N_Y_sel));
%         plot_single_result_panel(N_Y_sel, dtheta_vec, anly_sel, anly_field, ...
%             z_sem, cmap_adj, fill_alpha, pltsty, ttl, ylbl,  ...
%             lgdn_opts_on{(splt_cnt==nrows*ncols) + 1}{:})
        set(gca, 'tag', anly_field);
        
    end
end

despline('all');

linkaxes(findall(gcf,'type','axes'),'xy');

%%

figure; 

splt_cnt = 1;
prog_field = 'best'; 

for i_nx = 1:num_N_X
    N_X_sel = N_X_vec(i_nx); 
    sel_ind_N_X = find(combination_table.N_X == N_X_sel); 
    N_Y_sel = combination_table.N_Y(sel_ind_N_X); 
    
    [N_Y_sel, sorted_ind] = sort(N_Y_sel);
    
    anly_sel = analyses_aggr.(prog_field)(sel_ind_N_X);
    anly_sel = anly_sel(sorted_ind);
        
    for i_ny = 1:num_N_Y
        N_Y_val = N_Y_sel(i_ny);
        
        subplot(num_N_X, num_N_Y, splt_cnt); hold on;
        splt_cnt = splt_cnt + 1; 
        
        anly_sel_mat = anly_sel{i_ny}.prob.mean;
        
        image_with_strict_limits(log10(anly_sel_mat)'); 
%         daspect([1,2,1]);
        
%         colorbar; 
%         title(sprintf('%d X, %d Y', N_X_sel, N_Y_val))
        
        caxis([-10, 0]);
        if i_ny == 1
            ylabel({sprintf('%d X', N_X_sel), '$\Delta \theta$'});
        end
        if i_nx == num_N_X
            xlabel({'pattern', sprintf('%d Y', N_Y_val)});
        end
        
    end
end

%%

all_results_agg = cell(num_comb,1);
cnt_res = 1;
for i_nx = 1:num_N_X
    N_X_sel = N_X_vec(i_nx); 
    sel_ind_N_X = find(combination_table.N_X == N_X_sel); 
    N_Y_sel = combination_table.N_Y(sel_ind_N_X); 
    
    data_files = arrayfun(@(x) fullfile(data_path, sprintf('sim_%03d.mat', x)), sel_ind_N_X, 'uni', 0);
    
N_Y_vec = zeros(num_N_Y,1);

sm_wins = struct('entropy', 5, 'num_unq', 10);
net_fields = {'theta', 'W'}';

tic
for i = 1:num_N_Y
    dat = load(data_files{i}); 
    res = struct();
    
    N_Y_vec(i) = dat.configurations.N_Y;
    dtheta_vec = dat.configurations.dtheta_vec;
    
    L1_fields = fieldnames(dat.simulations{1}.net_config);
%     L2_fields = fieldnames(dat.simulations{1}.apply_dtheta.(L1_fields{1}));
    L2_fields = {'entropy', 'num_unq'}';
    
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
    
    all_results_agg{cnt_res} = res;
    cnt_res = cnt_res + 1;
end

toc


end

% all_results_agg = structarray_to_struct(vertcat(all_results_agg{:}));
% all_results_agg = structfun(@(S) structarray_to_struct(S, 0), all_results_agg, 'uni', 0);
% 
% 
% %%
%     L2_fields = {'entropy', 'num_unq'}';
% 
% 
% prog_fields = L1_fields;
% anly_fields = [{'net'}; L2_fields];
% 
% for i1 = 1:length(prog_fields)
%     F1 = prog_fields{i1};
% 
%     for i2 = 1:length(anly_fields)
%         F2 = anly_fields{i2};
%         all_results_agg.(F1).(F2) = structfun(@(x) horzcat(x{:})', ...
%             structarray_to_struct([all_results_agg.(F1).(F2){:}], 0), 'uni',0);
% 
%     end
% end

%%
figure; 