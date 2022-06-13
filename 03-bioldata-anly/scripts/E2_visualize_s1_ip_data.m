clc; clear; close all;
run start_up.m

%% Define paths
main_data_path = 'data/ip-data';

fig_path = 'figures/ip-data/sum'; 
if ~exist(fig_path, 'dir')
    mkdir(fig_path); 
end

PROCESS_path = main_data_path; 

%% Load data
load(fullfile(PROCESS_path, 'analysis-summary.mat'));

analysis_results = analysis_table.analysis;
representative_traces = analysis_table.representative;

cell_fullIDs = arrayfun(@(x) x.full_ID, analysis_table.info, 'uni', 0); 
exp_groups = arrayfun(@(x) x.Group, analysis_table.info, 'uni', 0);

%% Expgroups colors and reassignment 
unq_expgroups = {'cholinergic', 'electric', 'cholinergic paired'}; % customized order for colors instead of using `unique` 
expgroup_colors = flipud(return_colorbrewer('Set1', length(unq_expgroups))*0.95);
groupcolor_map = containers.Map(unq_expgroups, mat2colcell(expgroup_colors));

%% Load selection 
cell_selection = readtable(fullfile(PROCESS_path, 'cell-selection.csv'), 'PreserveVariableNames', true); 
cell_selection.base = cellfun(@(x) str2num(x), cell_selection.base, 'uni', 0);  %#ok<ST2NM>
cell_selection.post = cellfun(@(x) str2num(x), cell_selection.post, 'uni', 0);  %#ok<ST2NM>
selected_cellids = cell_selection.cell_id(cell_selection.selected == 1); 

select_conditions = contains(cell_fullIDs, selected_cellids);
select_analysis = analysis_results(select_conditions); 
select_cellids = cell_fullIDs(select_conditions);
select_num_cells = length(select_analysis);

%% Select measures
select_measures = {'Vthres_first', 'Vrest_1', 'dVthres_rest', 'num_spikes'}; 

meantime_fun = @(x,y,tr) mean(y(x >= min(tr) & x <= max(tr)), 'omitnan');

pooled_analysis = cell(select_num_cells, 1); 
pooled_expgroups = cell(select_num_cells, 1);

for i = 1:select_num_cells
    sel_obj = select_analysis(i);
    sel_id = select_cellids(i);
    t_vec = sel_obj.time_vec; 
    t_base = cell_selection.base{strcmp(cell_selection.cell_id, sel_id)};
    t_post = cell_selection.post{strcmp(cell_selection.cell_id, sel_id)};
    
    sel_obj.dVthres_rest = sel_obj.Vthres_first - sel_obj.Vrest_1;

    tmp_struct = struct; 
    for j = 1:length(select_measures)
        measure_j = select_measures{j}; 
        vec_j = sel_obj.(measure_j); 
        base_j = meantime_fun(t_vec, vec_j, t_base);
        post_j = meantime_fun(t_vec, vec_j, t_post);
        
        tmp_struct.([measure_j '_base']) = base_j; 
        tmp_struct.([measure_j '_post']) = post_j; 
        tmp_struct.([measure_j '_change']) = post_j - base_j;
        
        if any(strcmp(measure_j, {'dVthres_rest', 'num_spikes'}))
            tmp_struct.([measure_j '_percentchange']) = 100*(post_j - base_j) / base_j;
        end
    end
    pooled_analysis{i} = tmp_struct;
    pooled_expgroups{i} = exp_groups{strcmp(cell_fullIDs, sel_id)};
end

pooled_analysis = structarray_to_struct(vertcat(pooled_analysis{:})); 
fields_pooled = fieldnames(pooled_analysis);
values_pooled = struct2cell(pooled_analysis);

pooled_table = table(select_cellids, pooled_expgroups, 'VariableNames', {'CellID', 'ExpGroup'});
pooled_table = [pooled_table, struct2table(pooled_analysis)];

%% Latex names
latex_changed = struct(...
    'num_spikes', 'n_{spk}', ...
    'Vthres_first', 'V_T', ...
    'Vrest_1', 'V_R', ...
    'Rin', 'R_{\mathrm{in}}', ...
    'fAHP_Vm_mean', 'AHP', ...
    'dVthres_rest', '\Delta V_{TR}', ...
    'dVfAHP_rest', '\Delta V_{AHP}');
latex_suffix = struct(...
    'base', '0', ...
    'change', '\Delta', ...
    'percentchange',  '\%\Delta', ...
    'post', 'f'); 

latex_pooled = fields_pooled; 
latex_pooled_struct = structfun(@(x) '', pooled_analysis, 'uni', 0);
for i = 1:length(latex_pooled)
    raw_name = latex_pooled{i}; 
    sep_ind  = regexp(raw_name, '_'); sep_ind = sep_ind(end); 
    main_name = latex_changed.(raw_name(1:sep_ind-1)); 
    suff_name = latex_suffix.(raw_name(sep_ind+1:end)); 
    latex_pooled{i} = sprintf('$%s^{\\mathbf{%s}}$', main_name, suff_name);
    latex_pooled_struct.(raw_name) = latex_pooled{i};
end

%% Save for later use
save(fullfile(PROCESS_path, 'ip-select-data-analysis.mat'), ...
    'pooled_table', 'unq_expgroups', 'latex_pooled_struct');

%% Violin summary 
pairs_to_plt = {'Vrest_1_base', 'Vrest_1_post'; ...
                'Vthres_first_base', 'Vthres_first_post'; ...
                'Vrest_1_base', 'Vthres_first_base'; ...
                'Vrest_1_post', 'Vthres_first_post'; ...
                'dVthres_rest_base',  'dVthres_rest_post'};
               
ncols = size(pairs_to_plt,1);
scale_xaxis = 0.4; 

figure('units','normalized','position',[0,0,0.65,1]);

graphic_setdefault(35, ...
    'DefaultAxesMinorGridAlpha', 0.05, ...
    'DefaultAxesMinorGridLineStyle', '-', ...
    'DefaultTextInterpreter', 'latex', ...
    'DefaultLegendInterpreter', 'latex', ...
    'DefaultAxesTitleFontSize', 1.2, ...
    'DefaultAxesLabelFontSize', 1.2, ...
    'DefaultAxesLineWidth', 2, ...
    'DefaultLineLineWidth', 3);

ylim_1 = [-86, -29];
ylim_2 = [0, diff(ylim_1)];

for i = 1:ncols 
    subplot(1, ncols, i); hold on;
    field_1 = pairs_to_plt{i,1};
    field_2 = pairs_to_plt{i,2}; 
    
    xaxis_names = {latex_pooled_struct.(field_1), latex_pooled_struct.(field_2)};
    yv = [pooled_analysis.(field_1), pooled_analysis.(field_2)];
    xv = repmat([1,2], [size(yv,1),1]);
    xv = xv + scale_xaxis*(rand(size(xv)) - 0.5);
    vp = violinplot(yv, xaxis_names, ...
        'ShowData', false, 'ViolinColor', [0.6,0.6,0.6], ...
        'ShowNotches', false, 'BoxColor', [0.6,0.6,0.6], ...
        'ShowMean', true, 'Bandwidth', 3, ...
        'EdgeColor', [1,1,1]);
    
    for j = 1:length(unq_expgroups)
        expgroup_lgnd = unq_expgroups{j};
        expgroup_inds = strcmp(pooled_expgroups, expgroup_lgnd);
        xvj = xv(expgroup_inds,:);
        yvj = yv(expgroup_inds,:);
        clrj = expgroup_colors(j,:);
        plot(xvj',yvj', 'color', [clrj,0.25], 'linewidth', 2);
        scatter(xvj(:), yvj(:), 180, clrj, 'filled',...
            'markeredgecolor', 'none', ...
            'markerfacealpha',0.3, ...
            'displayname', expgroup_lgnd, ...
            'tag', 'showlegend');
    end
    
    set(gca, 'xtick', [1,2], ...
        'xticklabel', xaxis_names, ...
        'ticklabelinterpreter', 'latex');
    
    if i == 1 || i == ncols
        ylabel('$V_m$ (mV)');
    else
        set(gca, 'ycolor', 'none');
    end
    
    if i < ncols 
        set(gca, 'tag', 'linked');
        pos = get(gca, 'position'); 
        set(gca, 'position', pos + [-0.02,0,0,0]);
    else 
        ylim(ylim_2);
        pos = get(gca, 'position'); 
        set(gca, 'position', pos + [0.03,0,0.05,0]);
    end
    if i == 3
        [lgnd, icons, ~,~] = legend(findobj(gca, 'tag', 'showlegend'), 'Fontsize', 30);
        set(lgnd, 'Box', 'on', 'Color', [1,1,1]*0.97, 'Location', 'southwest');
        lgnd.Position = lgnd.Position + [-0.2, 0.01, 0, 0];
        arrayfun(@(x) set(x.Children, 'MarkerSize', 15), findall(icons, 'type', 'hggroup'));      
        
    end

end

linked_axes = findall(gcf, 'type', 'axes', 'tag', 'linked');
linkaxes(linked_axes, 'y'); 
despline('all'); 

ylim(linked_axes(1), ylim_1);

exportgraphics(gcf, fullfile(fig_path, 's1-violin.pdf'));
close; 

%% Plot change stats

pairs_to_plt = {...
    'dVthres_rest_change', 'num_spikes_change'; 
    'Vthres_first_change', 'num_spikes_change';
    'Vrest_1_change', 'num_spikes_change'};
xlbl_unit = ' (mV)'; 
ylbl_unit = ''; 

stripdollarsymbols = @(x) regexprep(x, '\$', '');
ncols = size(pairs_to_plt,1); 
nrows = 1; 
scale_expand_x = 0.2;
xlim_config_fun = @(xr) 5*round(xr/5);

graphic_setdefault(28, ...
    'DefaultAxesMinorGridAlpha', 0.05, ...
    'DefaultAxesMinorGridLineStyle', '-', ...
    'DefaultTextInterpreter', 'latex', ...
    'DefaultLegendInterpreter', 'latex', ...
    'DefaultAxesTitleFontSize', 1.2, ...
    'DefaultAxesLabelFontSize', 1.2, ...
    'DefaultAxesLineWidth', 2, ...
    'DefaultLineLineWidth', 3);


figure('units','normalized','position',[0,0,0.95,0.65]);

for i = 1:size(pairs_to_plt,1)
    subplot(nrows, ncols, i); hold on;
    
    ax = gca;
    ax.Position = ax.Position .* [1,1,1.25,0.8] + [-0.06 + (i-1)*0.04,0.06,0,0];
    
    field_1 = pairs_to_plt{i,1};
    field_2 = pairs_to_plt{i,2}; 
    
    axlbl_names = {...
        sprintf('%s%s', latex_pooled_struct.(field_1), xlbl_unit), ...
        sprintf('%s%s', latex_pooled_struct.(field_2), ylbl_unit)};
    xv = pooled_analysis.(field_1);
    yv = pooled_analysis.(field_2);
    
    % x-vector for fitting 
    range_xv = [min(xv);max(xv)];
    xv_f = [min(xv);max(xv)] + [-1;1]*scale_expand_x*diff(range_xv);
    
    % fitting whole
    mdl = fitlm(xv, yv);
    R2 = mdl.Rsquared.Ordinary;
    pval = mdl.coefTest;
    yv_f = mdl.predict(xv_f); 
    pval_txt = [regexprep(sprintf('%.2e', pval),'e', '\\times10^{'), '}'];
    ttl_txt = {...
        sprintf('$\\mathbf{%s \\sim %s + 1}$', ...
            stripdollarsymbols(latex_pooled_struct.(field_2)), ...
            stripdollarsymbols(latex_pooled_struct.(field_1))), ...
        sprintf('\\fontsize{%f}{0}\\selectfont $p = %s, R^2 = %.2f$', 8, pval_txt, R2)};
    
    plot(xv_f, yv_f, ...
        'linestyle', '-.', ...
        'color', [0.5,0.5,0.5], ...
        'linewidth', 5, ...
        'displayname', 'fit all', ...
        'tag', 'showlegend');
    
    % plot and fit line each group 
    for j = 1:length(unq_expgroups)
        expgroup_lgnd = unq_expgroups{j}; 
        expgroup_inds = strcmp(pooled_expgroups, expgroup_lgnd);
        xvj = xv(expgroup_inds,:);
        yvj = yv(expgroup_inds,:);
        clrj = expgroup_colors(j,:);
        scatter(xvj(:), yvj(:), 180, clrj, 'filled',...
            'markeredgecolor', 'none', ...
            'markerfacealpha',0.4, ...
            'displayname', expgroup_lgnd);
        
        lin_fit = polyfit(xvj,yvj,1);
        yvjf = lin_fit(1)*xv_f+lin_fit(2);
        plot(xv_f, yvjf, '-', 'color', [clrj, 0.2], 'linewidth', 5);
    end
    
    xline(0, 'k:', 'linewidth', 2);
    yline(0, 'k:', 'linewidth', 2);
    
    xlabel(axlbl_names{1});
    if i == 1
        ylabel(axlbl_names{2});
    end
    title(ttl_txt)
    xlim(xlim_config_fun(xv_f));
    
    if i == length(pairs_to_plt)
        [lgnd, icons, ~,~] = legend(findobj(gca, 'tag', 'showlegend'),'fontsize', 22);
        set(lgnd, 'Box', 'on', 'Color', [1,1,1]*0.97, 'Location', 'southeast');
        arrayfun(@(x) set(x.Children, 'MarkerSize', 15), findall(icons, 'type', 'hggroup'));        
    end

end

linkaxes(findall(gcf, 'type', 'axes'), 'y'); 
despline('all');
ylim([-6,10]);

annotation('textbox', 'units', 'normalized', 'position', [0.08,lgnd.Position(2),0.14,0.14], ...
    'string', ...
    {'$x^{\mathbf{\Delta}} = x^{\mathrm{f}} - x^{\mathrm{0}}$', ...
    ' (final $-$ initial)'}, ...
    'Interpreter', 'latex', 'fontsize', 25, 'LineWidth', 2, 'BackgroundColor', [1,1,1]*0.97,  ...
    'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left')

exportgraphics(gcf, fullfile(fig_path, 's1-change.pdf'));
close;

%% Sorted change based on dV_TR
[sorted_dVchange, sorted_ind] = sort(pooled_table.dVthres_rest_change);
sorted_rel_dVchange = pooled_table.dVthres_rest_percentchange(sorted_ind);

sorted_nspkchange = pooled_table.num_spikes_change(sorted_ind); 
sorted_rel_nspkchange = pooled_table.num_spikes_percentchange(sorted_ind); 

sorted_expgroup = pooled_table.ExpGroup(sorted_ind);

pseudo_xvec = 1:select_num_cells;

text_dVTR_ypos = -sign(sorted_dVchange) * 6;
combined_dVTR_text = arrayfun(@(x,y) sprintf('%.2f mV, \\color[rgb]{0.6,0.6,0.6}%.1f %%', x,y), ...
    sorted_dVchange, sorted_rel_dVchange, 'uni', 0);

text_nspk_ypos = -sign(sorted_nspkchange) * 5;
combined_nspk_text = arrayfun(@(x,y) sprintf('%.2f, \\color[rgb]{0.6,0.6,0.6}%.1f %%', x,y), ...
    sorted_nspkchange, sorted_rel_nspkchange, 'uni', 0);

stem_args = {'linewidth', 7, 'linestyle', '-', 'marker', 'none'};

figure('units','normalized','position',[0,0,1,1]); 

subplot(3,1,[1,2]);
hold on; 

stem(sorted_dVchange,'k', stem_args{:});
ylim([-12, 18]);
ylabel(sprintf('%s (mV)', latex_pooled_struct.dVthres_rest_change));

text(pseudo_xvec, text_dVTR_ypos, combined_dVTR_text, ...
    'Interpreter','tex', 'HorizontalAlignment', 'center',...
    'FontSize', 15, 'Rotation', 90)

for i = 1:length(unq_expgroups)
    exp_group_ith = unq_expgroups{i};
    loc2plt = strcmp(sorted_expgroup, exp_group_ith);
    
    scatter(pseudo_xvec(loc2plt), 17 * ones(1,sum(loc2plt)), ...
        250, expgroup_colors(i,:), ...
        'filled', 'MarkerEdgeColor', 'none', ...
        'displayname', exp_group_ith, 'tag', 'showlegend');
end

set(gca, 'ytick', -12:6:18);

[lgnd, icons, ~, ~] = legend(findobj(gca, 'tag', 'showlegend'),'fontsize', 19);
set(lgnd, 'Box', 'on', 'Color', [1,1,1]*0.97, 'Location', 'North');
lgnd.Position = lgnd.Position + [0.2,-0.05,0,0]; 
arrayfun(@(x) set(x.Children, 'MarkerSize', 15), findall(icons, 'type', 'hggroup'));

yyaxis right;
stem(pseudo_xvec + 0.2, sorted_rel_dVchange, 'color', [0.6,0.6,0.6], stem_args{:});
ylim([-100, 150]);
ylabel(sprintf('%s \\%%', latex_pooled_struct.dVthres_rest_percentchange));
set(gca, 'ycolor', [0.6,0.6,0.6]);

xlim([0,select_num_cells+1]);
hide_only_axis('x');
title(['S1 cells sorted by ' latex_pooled_struct.dVthres_rest_change]); 

subplot(3,1,3);
hold on; 

stem(sorted_nspkchange,'k', stem_args{:});
ylabel(sprintf('%s', latex_pooled_struct.num_spikes_change));
text(pseudo_xvec, text_nspk_ypos, combined_nspk_text, ...
    'Interpreter','tex', 'HorizontalAlignment', 'center',...
    'FontSize', 15, 'Rotation', 90)

text_ypos_cellids = 10 * ones(1,select_num_cells);
stem(pseudo_xvec + 0.1, text_ypos_cellids, ':', 'color', [0.3,0.1,0.6,0.5]);
sorted_cellids = pooled_table.CellID(sorted_ind); 
text(pseudo_xvec, text_ypos_cellids, sorted_cellids, ...
    'Interpreter','none', 'FontAngle', 'italic', ...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', ...
    'FontSize', 8, 'Rotation', 60, ...
    'BackgroundColor', [0.4,0.2,0.7], 'Color', 'w')

ylim([-5, 10]);
yyaxis right;
stem(pseudo_xvec + 0.2, sorted_rel_nspkchange, 'color', [0.6,0.6,0.6], stem_args{:});

set(gca, 'ycolor', [0.6,0.6,0.6]);
ylabel(sprintf('%s \\%%', latex_pooled_struct.num_spikes_percentchange));
ylim([-100, 200]);

xlim([0,select_num_cells+1]);
hide_only_axis('x');

linkaxes(findall(gcf, 'type', 'axes'), 'x');
exportgraphics(gcf, fullfile(fig_path, 's1-sorteddVchange.pdf'));
close;

%% Sorted change based on V_T
[sorted_VTchange, sorted_ind] = sort(pooled_table.Vthres_first_change);
sorted_dVchange = pooled_table.dVthres_rest_change(sorted_ind);

sorted_nspkchange = pooled_table.num_spikes_change(sorted_ind); 
sorted_rel_nspkchange = pooled_table.num_spikes_percentchange(sorted_ind); 

sorted_expgroup = pooled_table.ExpGroup(sorted_ind);

pseudo_xvec = 1:select_num_cells;

text_V_ypos = -sign(sorted_VTchange) * 5;
combined_V_text = arrayfun(@(x,y) sprintf('%.2f mV, \\color[rgb]{0.6,0.6,0.6}%.2f mV', x,y), ...
    sorted_VTchange, sorted_dVchange, 'uni', 0);

text_nspk_ypos = -sign(sorted_nspkchange) * 5;
combined_nspk_text = arrayfun(@(x,y) sprintf('%.2f, \\color[rgb]{0.6,0.6,0.6}%.1f %%', x,y), ...
    sorted_nspkchange, sorted_rel_nspkchange, 'uni', 0);

stem_args = {'linewidth', 7, 'linestyle', '-', 'marker', 'none'};

figure('units','normalized','position',[0,0,1,1]); 

subplot(3,1,[1,2]);
hold on; 

stem(sorted_VTchange,'k', stem_args{:});
ylim([-12, 12]);
ylabel(sprintf('%s (mV)', latex_pooled_struct.Vthres_first_change));

text(pseudo_xvec-0.1, text_V_ypos, combined_V_text, ...
    'Interpreter','tex', 'HorizontalAlignment', 'center',...
    'FontSize', 15, 'Rotation', 90)

for i = 1:length(unq_expgroups)
    exp_group_ith = unq_expgroups{i};
    loc2plt = strcmp(sorted_expgroup, exp_group_ith);
    
    scatter(pseudo_xvec(loc2plt), 12 * ones(1,sum(loc2plt)), ...
        250, expgroup_colors(i,:), ...
        'filled', 'MarkerEdgeColor', 'none', ...
        'displayname', exp_group_ith, 'tag', 'showlegend');
end

set(gca, 'ytick', -12:6:18);

[lgnd, icons, ~, ~] = legend(findobj(gca, 'tag', 'showlegend'),'fontsize', 19);
set(lgnd, 'Box', 'on', 'Color', [1,1,1]*0.97, 'Location', 'North');
lgnd.Position = lgnd.Position + [0.2,-0.05,0,0]; 
arrayfun(@(x) set(x.Children, 'MarkerSize', 15), findall(icons, 'type', 'hggroup'));

yyaxis right;
stem(pseudo_xvec + 0.2, sorted_dVchange, 'color', [0.6,0.6,0.6], stem_args{:});
ylim([-15,15]);
ylabel(sprintf('%s (mV)', latex_pooled_struct.dVthres_rest_change));
set(gca, 'ycolor', [0.6,0.6,0.6]);

xlim([0,select_num_cells+1]);
hide_only_axis('x');
title(['S1 cells sorted by ' latex_pooled_struct.Vthres_first_change]); 

subplot(3,1,3);
hold on; 

stem(sorted_nspkchange,'k', stem_args{:});
ylabel(sprintf('%s', latex_pooled_struct.num_spikes_change));
text(pseudo_xvec, text_nspk_ypos, combined_nspk_text, ...
    'Interpreter','tex', 'HorizontalAlignment', 'center',...
    'FontSize', 15, 'Rotation', 90)

text_ypos_cellids = 10 * ones(1,select_num_cells);
stem(pseudo_xvec + 0.1, text_ypos_cellids, ':', 'color', [0.3,0.1,0.6,0.5]);
sorted_cellids = pooled_table.CellID(sorted_ind); 
text(pseudo_xvec, text_ypos_cellids, sorted_cellids, ...
    'Interpreter','none', 'FontAngle', 'italic', ...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', ...
    'FontSize', 8, 'Rotation', 60, ...
    'BackgroundColor', [0.4,0.2,0.7], 'Color', 'w')

ylim([-5, 10]);
yyaxis right;
stem(pseudo_xvec + 0.2, sorted_rel_nspkchange, 'color', [0.6,0.6,0.6], stem_args{:});

set(gca, 'ycolor', [0.6,0.6,0.6]);
ylabel(sprintf('%s \\%%', latex_pooled_struct.num_spikes_percentchange));
ylim([-100, 200]);

xlim([0,select_num_cells+1]);
hide_only_axis('x');

linkaxes(findall(gcf, 'type', 'axes'), 'x');
exportgraphics(gcf, fullfile(fig_path, 's1-sortedVTchange.pdf'));
close;

%% Chi2 test of independence of dV_TR

cutoff_cmap = return_colorbrewer('Dark2', 3);
conf_levels = [0.1,0.01,0.05]; 

chi2plt_opts = {...
    struct(...
        'var_name', 'dVthres_rest_change', ...
        'var_cutoff_vec', linspace(0.01,2,200), ...
        'var_unit', 'mV'), ...
    struct(...
        'var_name', 'Vthres_first_change', ...
        'var_cutoff_vec', linspace(0.01,3,200), ...
        'var_unit', 'mV'), ...
    struct(...
        'var_name', 'dVthres_rest_percentchange', ...
        'var_cutoff_vec', linspace(0.05,6,200), ...
        'var_unit', '\%'), ...
    struct(...
        'var_name', 'Vrest_1_change', ...
        'var_cutoff_vec', linspace(0.01,2,200), ...
        'var_unit', 'mV')
};

cutoff_types = [-1,0,1];

nrow = 2;
ncol = 2;

graphic_setdefault(18, ...
    'DefaultAxesMinorGridAlpha', 0.05, ...
    'DefaultAxesMinorGridLineStyle', '-', ...
    'DefaultTextInterpreter', 'latex', ...
    'DefaultLegendInterpreter', 'latex', ...
    'DefaultAxesTitleFontSize', 1.4, ...
    'DefaultAxesLabelFontSize', 1.1, ...
    'DefaultAxesLineWidth', 2, ...
    'DefaultLineLineWidth', 3);

figure; 
for i = 1:length(chi2plt_opts)
    subplot(nrow,ncol,i); hold on;
    chi2_opt = chi2plt_opts{i};
    var_name = chi2_opt.var_name;
    latex_name = latex_pooled_struct.(var_name);
    var_cutoffs = chi2_opt.var_cutoff_vec; 
    var_unit = chi2_opt.var_unit;
    
    for j = 1:length(cutoff_types)
        cutoff_type = cutoff_types(j);
        
        p_vals = arrayfun(@(var_cutoff) chi2test_groupeffect(pooled_table,var_name, var_cutoff, cutoff_type), var_cutoffs);
        plot(var_cutoffs, p_vals,'color', cutoff_cmap(j,:), ...
            'displayname', cutofftype2signgroup(cutoff_type, true), ...
            'tag', 'showlegend');
    end
    
    arrayfun(@(y) yline(y, ':k', 'linewidth', 2), conf_levels);
    xlabel(sprintf('cutoff value (%s)', var_unit));
    ylabel('$p$ value');
    
    if i == 1
        arrayfun(@(y) text(max(var_cutoffs)*0.9, y, sprintf('\\textit{%.2f}', y), ...
            'FontSize', 20, 'VerticalAlignment', 'top'), conf_levels);
    end
    
    if i == 2
        lgnd = legend(findobj(gca, 'tag', 'showlegend'), 'fontsize', 12);
        title(lgnd, 'Sign groups');
        set(lgnd, 'Box', 'on');
    end
    
    title(sprintf('split by %s', latex_name));
    ylim([2e-3,1]);
    ax_pos = get(gca, 'position');
    ax_pos = ax_pos.*[1,1,1,0.9] + [0,-0.02,0,0];
    set(gca, 'yscale', 'log', 'position', ax_pos);

end

annotation('textbox', 'units', 'normalized', 'position', [0,0.95,1.0,0.05], ...
    'string', '$\chi^2$ test of independence for [\textit{sign-group} $\times$ \textit{exp-group}]', ...
    'Interpreter', 'latex', 'fontsize', 30, 'Linestyle', 'none', ...
    'VerticalAlignment', 'top', 'HorizontalAlignment', 'center')

pause(2);
savefig(gcf, fullfile(fig_path, 's1-chi2-Vm'));
exportgraphics(gcf, fullfile(fig_path, 's1-chi2-Vm.pdf'));
close;


%% Demo for split group 
signgroup_cmap = [0.05,0.05,0.05;0.4,0.4,0.4;0.8,0.8,0.8];
txtgroup_cmap  = [0.0,0.0,0.0;0.98,0.98,0.98;0.95,0.95,0.95];

% demo_chi2_opts = {...
%     struct(...
%         'var_name', 'dVthres_rest_change', ...
%         'var_cutoff', dVTR_change_cutoff, ...
%         'var_additional_null', @(T) T.num_spikes_change < 0 & T.dVthres_rest_change < 0, ...
%         'cutoff_types', -1, ...
%         'var_unit', 'mV');
% };

demo_chi2_opts = {...
    struct(...
        'var_name', 'dVthres_rest_percentchange', ...
        'var_cutoff', 1, ...               
        'var_additional_null', @(T) T.num_spikes_change < 0 & T.dVthres_rest_change < 0, ...
        'cutoff_types', -1, ...
        'var_unit', '%');
};

fig_file_name = 's1-chi2-pie.pdf'; 

graphic_setdefault(20, ...
    'DefaultAxesMinorGridAlpha', 0.05, ...
    'DefaultAxesMinorGridLineStyle', '-', ...
    'DefaultTextInterpreter', 'latex', ...
    'DefaultLegendInterpreter', 'latex', ...
    'DefaultAxesTitleFontSize', 1.4, ...
    'DefaultAxesLabelFontSize', 1.1, ...
    'DefaultAxesLineWidth', 2, ...
    'DefaultLineLineWidth', 3);

for demo_chi2_opt = demo_chi2_opts
    demo_chi2_opt = demo_chi2_opt{:}; %#ok<FXSET>
    var_name = demo_chi2_opt.var_name;
    var_cutoff = demo_chi2_opt.var_cutoff;
    var_unit = demo_chi2_opt.var_unit; 
    var_latex = latex_pooled_struct.(var_name);
    
    if isfield(demo_chi2_opt, 'var_additional_null')
        additional_args = {demo_chi2_opt.var_additional_null}; 
    end 
    
    for cutoff_type = demo_chi2_opt.cutoff_types
        
        fig_name = sprintf('%s-[var=%s]-[cutoff=%.2f]-[type=%d].pdf', fig_file_name, var_name, var_cutoff, cutoff_type);
        
        [pval, chi2aux] = chi2test_groupeffect(pooled_table, var_name, var_cutoff, cutoff_type, additional_args{:});
        
        sign_group_data = chi2aux.data;
        pie_data = sign_group_data.Variables;
        pie_titles = sign_group_data.Properties.RowNames;
        pie_groups = sign_group_data.Properties.VariableNames;      
        
        pie_group_lbls = cellfun(@(x) sprintf('%s threshold shift', x) , pie_groups, 'uni', 0);
        figure('units', 'normalized', 'position', [0,0,0.95,0.4]);
        
        tiledlayout(1,length(pie_titles), 'TileSpacing','none');
        
        for i = 1:length(pie_titles)
            ax = nexttile;
            data_ith = pie_data(i,:);
            n_sample_ith = sum(pie_data(i,:));
            pie(ax, data_ith, '$\\mathbf{%.0f\\%%}$');
            
            txt_objs = findobj(ax,'type','text');
            for j = 1:length(txt_objs)
                ptxt = txt_objs(j); 
                txt_color = txtgroup_cmap(j,:);   
                set(ptxt, 'FontSize', 22, 'Color', txt_color, ...
                    'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
                ptxt.Position = 0.65 * ptxt.Position;
            end
            
            set(findobj(ax,'type','patch'), 'LineWidth', 4, 'EdgeColor', 'w');
            ax.Colormap = signgroup_cmap;
            ttl = title(ax,sprintf('%s (n=%d)', pie_titles{i}, n_sample_ith));
            set(ttl,'unit', 'normalized');
            ttl.Position(2) = ttl.Position(2) - 1.1;
        end
        
        lgnd = legend(pie_group_lbls, 'Fontsize', 24);
        lgnd.Layout.Tile = 'east';
        
        lgnd_ttl = sprintf('$\\chi^2$ test $p$-value = %.3g', pval);
        % lgnd_ttl = {lgnd_ttl, sprintf('split by %s $\\sim %.2f$', var_latex, var_cutoff)};
        title(lgnd, lgnd_ttl, 'Interpreter', 'latex');
        
%         exportgraphics(gcf, fullfile(fig_path, fig_name));
%         close;
    end
end

%% Helper functions 
function [pval,aux] = chi2test_groupeffect(tbl, var_name, var_cutoff, cutoff_type, var_additional_null)
    
    SG_add_null = zeros(height(tbl),1);
    if exist('var_additional_null', 'var')
        SG_add_null = var_additional_null(tbl);
    end
    
    [EG, EG_cat] = grp2idx(tbl.ExpGroup);
    X = tbl.(var_name);
    SG = sign(X);

    var_cutoff = abs(var_cutoff);

    switch cutoff_type
        case 0
            SG(abs(X) <= var_cutoff | SG_add_null) = 0;
        case 1
            SG(X <= var_cutoff | SG_add_null) = -1;
        case -1
            SG(X >= -var_cutoff | SG_add_null) = 1;
        otherwise
            error('Only -1,0,1 are accepted');
    end
    sign_group = cutofftype2signgroup(cutoff_type);
    
    [SG, ~] = grp2idx(SG);

    P_X = histcounts2(EG, SG, 'BinMethod', 'integers');

    [~,pval,chi2,df] = chi2tod(P_X);

    if nargout == 1
        return
    end

    aux = struct('data', struct, 'stat', struct);
    aux.data = array2table(P_X, 'VariableNames', sign_group, 'RowNames', EG_cat);
    aux.stat = struct('pval', pval, 'chi2', chi2', 'df', df);

end

function res = cutofftype2signgroup(cutoff_type, is_legend)
    if ~exist('is_legend', 'var'), is_legend = false; end

    switch cutoff_type
        case 0
            res = {'negative','0','positive'};            
        case 1
            res = {'non-positive','positive'};
        case -1
            res = {'negative','non-negative'};
        otherwise
            error('Only -1,0,1 are accepted');
    end
    
    if ~is_legend, return; end
    res = sprintf('(%s)', strjoin(res, ', '));
    
end

function res = append_varname_signgroup(sign_group, var_name, var_cutoff, var_unit, val_fmt)
    if ~exist('var_unit', 'var')
        var_unit = ''; 
    else 
        var_unit = sprintf(' %s', var_unit);
    end
    
    if ~exist('val_fmt', 'var'), val_fmt = '%.1f'; end
    var_name = regexprep(var_name, '\$', ''); 
    
    switch sign_group
        case '0'
            var_fmt = '\\left| %s \\right| \\leq';
        case 'negative'
            var_fmt = '%s < -';
        case 'positive'
            var_fmt = '%s >';
        case 'non-positive'
            var_fmt = '%s \\leq';
        case 'non-negative'
            var_fmt = '%s \\geq -';
    end
    fmt = ['$', var_fmt, val_fmt, '%s$'];
    res = sprintf(fmt, var_name, var_cutoff, var_unit);

end