clc; clear; close all;
run start_up.m

graphic_setdefault(25, ...
    'DefaultAxesMinorGridAlpha', 0.05, ...
    'DefaultAxesMinorGridLineStyle', '-', ...
    'DefaultTextInterpreter', 'latex', ...
    'DefaultLegendInterpreter', 'latex', ...
    'DefaultAxesTitleFontSize', 1.1, ...
    'DefaultAxesLabelFontSize', 1.1, ...
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
select_measures = {'Vthres_first', 'Vrest_1', 'dVthres_rest'}; 

basetime_fun = @(x,y,tb) mean(y(x >= min(tb) & x <= max(tb)), 'omitnan');
posttime_fun = @(x,y,tp) mean(y(x >= min(tp) & x <= max(tp)), 'omitnan');

pooled_analysis = cell(select_num_cells, 1); 
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
        tmp_struct.([measure_j '_base']) = base_j; 
        tmp_struct.([measure_j '_post']) = post_j; 
        tmp_struct.([measure_j '_change']) = post_j - base_j;
    end
    pooled_analysis{i} = tmp_struct;
end

pooled_analysis = structarray_to_struct(vertcat(pooled_analysis{:})); 
fields_pooled = fieldnames(pooled_analysis);
values_pooled = struct2cell(pooled_analysis);

%% Latex names
latex_changed = struct(...
    'num_spikes', 'n_{spk}', ...
    'Vthres_first', 'V_T', ...
    'Vrest_1', 'V_R', ...
    'Rin', 'R_i', ...
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

%% Violin summary 
pairs_to_plt = {'Vrest_1_base', 'Vrest_1_post'; ...
                'Vthres_first_base', 'Vthres_first_post'; ...
                'Vrest_1_base', 'Vthres_first_base'; ...
                'Vrest_1_post', 'Vthres_first_post'; ...
                'dVthres_rest_base',  'dVthres_rest_post'};
               
ncols = size(pairs_to_plt,1);
scale_xaxis = 0.2; 

figure('units','normalized','position',[0,0,0.6,0.7]);
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
        'ShowNotches', true, 'BoxColor', [0.6,0.6,0.6], ...
        'ShowMean', true, 'Bandwidth', 3, ...
        'EdgeColor', [1,1,1]);
    
    plot(xv',yv', 'color', [0.1,0.1,0.1,0.15], 'linewidth', 0.75);
    scatter(xv(:), yv(:), 25, [0.1,0.1,0.1], 'filled',...
        'markeredgecolor', 'none', 'markerfacealpha',0.3);
    set(gca, 'xtick', [1,2], 'xticklabel', xaxis_names, 'ticklabelinterpreter', 'latex');
    
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
        pos = get(gca, 'position'); 
        set(gca, 'position', pos + [0.03,0,0,0]);
    end

end

linkaxes(findall(gcf, 'type', 'axes', 'tag', 'linked'), 'y'); 
despline('all'); 

% export_fig(fullfile(fig_path, 'dan-violin'), '-r200', '-p0.02')

%% Plot scatters 
pairs_to_plt = {'Vrest_1', 'Vthres_first'; ...
                'Vrest_1', 'dVthres_rest'};
unit_str = '(mV)';
ylim_cust = {[-55,-34],[-1, 35]};

cust_lims = struct('x', [-76,-59.5], 'y', {[-50.4767  -34.3240], [16.2592   33.0856]});

nrows = size(pairs_to_plt, 1); 
base_color = [0.1,0.1,0.1]; 
post_color = [0.8,0.2,0.2]; 

figure('units','normalized','position',[0,0,0.65,0.5]);
for i = 1:nrows
    pref_1 = pairs_to_plt{i,1}; 
    pref_2 = pairs_to_plt{i,2}; 
    
    x_base = pooled_analysis.([pref_1 '_base']);
    x_post = pooled_analysis.([pref_1 '_post']);
    x_name = sprintf('$%s$ %s', latex_changed.(pref_1), unit_str);
    
    y_base = pooled_analysis.([pref_2 '_base']);
    y_post = pooled_analysis.([pref_2 '_post']);
    y_name = sprintf('$%s$ %s', latex_changed.(pref_2), unit_str);
    
    lgdn_obj = gobjects(2,1);
    subplot(1,nrows,i); hold on; 
    plot([x_base x_post]', [y_base y_post]', '-', 'color', [0.1,0.1,0.1,0.2], 'linewidth', 2, 'marker', 'none')
    lgdn_obj(1) = scatter(x_base, y_base, 100, base_color, 'filled', 'MarkerEdgeColor', 'none', ...
        'MarkerFaceAlpha', 0.7,  'displayname', 'base');
    lgdn_obj(2) = scatter(x_post, y_post, 100, post_color, 'filled', 'MarkerEdgeColor', 'none', ...
        'MarkerFaceAlpha', 0.7,  'displayname', 'post');
    if i == 1
        plot(ylim_cust{i}, ylim_cust{i}, '--k', 'linewidth', 2);
    end
    
    cust_rect = cust_lims(i); 
    rectangle('Position', ...
        [cust_rect.x(1),cust_rect.y(1),diff(cust_rect.x),diff(cust_rect.y)], ...
        'LineWidth', 2, 'LineStyle', ':', 'EdgeColor', [0.7,0.7,0.7], 'tag', 'deleteable');
    daspect([1,1,1]);
    if i == nrows
        [~,lgnd_child] = legend(lgdn_obj, 'location', 'southwest');
        arrayfun(@(x) set(x.Children, 'markersize', 15), lgnd_child);
    end
    xlabel(x_name);
    ylim(ylim_cust{i});
    ylabel(y_name); 
end

linkaxes(findall(gcf, 'type', 'axes'), 'x'); 
xlim([-80,-40]);
despline('all'); 

% export_fig(fullfile(fig_path, 'dan-scatter'), '-r200', '-p0.02')

delete(findall(gcf, 'tag', 'deleteable')); 
for i = 1:nrows
    subplot(1,nrows,i); 
    xlim(cust_lims(i).x);
    ylim(cust_lims(i).y);
end
% export_fig(fullfile(fig_path, 'dan-scatter-zoomin'), '-r200', '-p0.02')

%% Conversion 
