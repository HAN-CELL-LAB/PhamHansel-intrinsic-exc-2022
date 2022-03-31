clc; clear; close all;
run start_up.m

graphic_setdefault(25, ...
    'DefaultAxesMinorGridAlpha', 0.05, ...
    'DefaultAxesMinorGridLineStyle', '-', ...
    'DefaultTextInterpreter', 'latex', ...
    'DefaultLegendInterpreter', 'latex', ...
    'DefaultAxesTitleFontSize', 1.15, ...
    'DefaultAxesLabelFontSize', 1.15, ...
    'DefaultAxesLineWidth', 2, ...
    'DefaultLineLineWidth', 3);

%% Define paths
main_data_path = 'data/gill-data/Layer II.III S1 BC IP';

fig_path = 'figures/gill-data/Layer II.III S1 BC IP/sum'; 
if ~exist(fig_path, 'dir')
    mkdir(fig_path); 
end

INFO_path = fullfile(main_data_path, 'info');
PROCESS_path = fullfile(main_data_path, 'processed');

%% Load data
load(fullfile(PROCESS_path, 'analysis-summary.mat'));

analysis_results = analysis_table.analysis;
representative_traces = analysis_table.representative;

cell_fullIDs = arrayfun(@(x) x.full_ID, analysis_table.info, 'uni', 0); 
exp_groups = arrayfun(@(x) x.Group, analysis_table.info, 'uni', 0);

%% Expgroups colors and reassignment 

expgroup_legends = table(...
    {'WT Somatic','WT Synaptic','WT oxo-m','WT oxo-m Somatic'}', ...
    {'elec-act','elec-act','cholingergic','chol + elec-act'}', ...
    'VariableNames', {'name','legend'});
unq_expgroups = {'cholingergic', 'elec-act', 'chol + elec-act'};

for i = 1:size(expgroup_legends,1)
    exp_groups(strcmp(expgroup_legends(i,:).name, exp_groups)) = expgroup_legends(i,:).legend;
end

expgroup_colors = return_colorbrewer('Set1', length(unq_expgroups))*0.95;
expgroup_colors = flipud(expgroup_colors);

%% Load selection 
cell_selection = readtable(fullfile(PROCESS_path, 'cell-selection.csv'), 'PreserveVariableNames', true); 
cell_selection.base = cellfun(@(x) str2num(x), cell_selection.base, 'uni', 0);  %#ok<ST2NM>
cell_selection.post = cellfun(@(x) str2num(x), cell_selection.post, 'uni', 0);  %#ok<ST2NM>
selected_cellids = cell_selection.cell_id(cell_selection.selected == 1); 

select_conditions = contains(cell_fullIDs, selected_cellids);

select_analysis = analysis_results(select_conditions);
% select_representative = representative_traces(select_conditions);
% select_cellidlatex = cellids_latex(select_conditions); 
select_cellids = cell_fullIDs(select_conditions);
% select_fignames = fignames(select_conditions); 

select_num_cells = length(select_analysis);

%% Select measures
select_measures = {'Vthres_first', 'Vrest_1', 'dVthres_rest', 'Rin', 'num_spikes'}; 

basetime_fun = @(x,y,tb) mean(y(x >= min(tb) & x <= max(tb)), 'omitnan');
posttime_fun = @(x,y,tp) mean(y(x >= min(tp) & x <= max(tp)), 'omitnan');

pooled_analysis = cell(select_num_cells, 1); 
pooled_expgroups = cell(select_num_cells, 1);

for i = 1:select_num_cells
    sel_obj = select_analysis(i);
    sel_id = select_cellids(i);
    t_vec = sel_obj.time_vec; 
    t_base = cell_selection.base{strcmp(cell_selection.cell_id, sel_id)};
    t_post = cell_selection.post{strcmp(cell_selection.cell_id, sel_id)};
    
    sel_obj.dVthres_rest = sel_obj.Vthres_first - sel_obj.Vrest_1;
    sel_obj.dVfAHP_rest = sel_obj.fAHP_Vm_mean - sel_obj.Vrest_1;
   
    tmp_struct = struct; 
    for j = 1:length(select_measures)
        measure_j = select_measures{j}; 
        vec_j = sel_obj.(measure_j); 
        base_j = basetime_fun(t_vec, vec_j, t_base);
        post_j = posttime_fun(t_vec, vec_j, t_post);
%         tmp_struct.([measure_j '_base']) = base_j; 
%         tmp_struct.([measure_j '_post']) = post_j; 
        tmp_struct.([measure_j '_change']) = post_j - base_j;
    end
    pooled_analysis{i} = tmp_struct;
    pooled_expgroups{i} = exp_groups{strcmp(cell_fullIDs, sel_id)};
end

pooled_analysis = structarray_to_struct(vertcat(pooled_analysis{:})); 
fields_pooled = fieldnames(pooled_analysis);
values_pooled = struct2cell(pooled_analysis);

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


%%
pairs_to_plt = {...
    'dVthres_rest_change', 'num_spikes_change'; 
    'Vthres_first_change', 'num_spikes_change';
    'Vrest_1_change', 'num_spikes_change'; };
%     'Rin_change', 'num_spikes_change'};

stripdollarsymbols = @(x) regexprep(x, '\$', '');
ncols = size(pairs_to_plt,1); 
nrows = 1; 
scale_expand_x = 0.2;

% figure('units','normalized','position',[0,0,0.8,1]);
figure('units','normalized','position',[0,0,1,0.6]);

for i = 1:size(pairs_to_plt,1)
    subplot(nrows, ncols, i); hold on;
    
    ax = gca;
    ax.Position = ax.Position .* [1,1,1.25,0.8] + [-0.07 + (i-1)*0.025,0.06,0,0];
    
    field_1 = pairs_to_plt{i,1};
    field_2 = pairs_to_plt{i,2}; 
    
    axlbl_names = {latex_pooled_struct.(field_1), latex_pooled_struct.(field_2)};
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
        sprintf('$\\mathbf{%s \\sim %s + 1}$', stripdollarsymbols(axlbl_names{1}), stripdollarsymbols(axlbl_names{2})), ...
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
            'displayname', expgroup_lgnd, ...
            'tag', 'showlegend');
        
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
    xlim(xv_f);
    if i == length(pairs_to_plt)
        [lgnd, icons, dsc,scd] = legend(findobj(gca, 'tag', 'showlegend'),'fontsize', 20);
        set(lgnd, 'Box', 'on', 'Color', [1,1,1]*0.97, 'Location', 'southwest');
        arrayfun(@(x) set(x.Children, 'MarkerSize', 15), findall(icons, 'type', 'hggroup'));        
    end

end

linkaxes(findall(gcf, 'type', 'axes'), 'y'); 
despline('all'); 
ylim([-6,10]);

annotation('textbox', 'units', 'normalized', 'position', [0.37,lgnd.Position(2),0.16,0.135], ...
    'string', ...
    {'$x^{\mathbf{\Delta}} = x^{\mathrm{f}} - x^{\mathrm{0}}$', ...
    'change = final - initial'}, ...
    'Interpreter', 'latex', 'fontsize', 20, 'LineWidth', 2, 'BackgroundColor', [1,1,1]*0.97,  ...
    'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left')

exportgraphics(gcf, fullfile(fig_path, 'test-dan-changestat-allexpgroup.pdf'))
% export_fig(fullfile(fig_path, 'dan-changestat-cholact'), '-r200', '-p0.02')
% 
% export_fig(fullfile(fig_path, 'dan-changestat-allexpgroup'), '-r200', '-p0.02', '-pdf')
% % export_fig(fullfile(fig_path, 'dan-violin-allexpgroup'), '-p0.02', '-svg')
% close; 

