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
    {'elec','elec','chol','chol-elec'}', ...
    'VariableNames', {'name','legend'});

unq_expgroups = {'chol', 'elec', 'chol-elec'};

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
select_cellids = cell_fullIDs(select_conditions);

select_num_cells = length(select_analysis);

%% Select measures
select_measures = {'Vthres_first', 'Vrest_1', 'dVthres_rest', 'num_spikes'}; 

renamed_measures = struct(...
    'Vthres_first', 'V_T', ...
    'Vrest_1', 'V_R', ...
    'dVthres_rest', 'dV_TR', ...
    'num_spikes', 'n_spk');

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
        
        measure_j = renamed_measures.(measure_j);
        tmp_struct.([measure_j '_base']) = base_j; 
        tmp_struct.([measure_j '_post']) = post_j; 
        tmp_struct.([measure_j '_change']) = post_j - base_j;
        
        if any(strcmp(measure_j, {'dV_TR', 'n_spk'}))
            tmp_struct.([measure_j '_percent_change']) = 100*(post_j - base_j) / base_j;
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

pooled_table.exp_elec = contains(pooled_table.ExpGroup,'elec');
pooled_table.exp_chol = contains(pooled_table.ExpGroup,'chol');

exp_group_count_table = groupcounts(pooled_table, 'ExpGroup');

pooled_vars = pooled_table.Properties.VariableNames;
sel_pooled_vars = pooled_vars(cellfun(@(x) ~isempty(x), regexp(pooled_vars, '(ExpGroup|.*_change)', 'match')));
changed_table = pooled_table(:,sel_pooled_vars);


%% 


y_vars = {'V_T_change', 'V_R_change', 'dV_change', 'n_spk_change'};
group_vars = {'none', 'elec', 'chol'};
res = struct; 
res_pval = struct;

for i = 1:length(y_vars)
    y_var = y_vars{i};
    res.(y_var) = struct; 
    res_pval.(y_var) = struct;
    for j = 1:length(group_vars)
        group_var = group_vars{j};
        res.(y_var).(group_var) = compare_expgroup_anov1(changed_table, 'ExpGroup', y_var, group_var);
        res_pval.(y_var).(group_var) = res.(y_var).(group_var).p_val_anova;
    end
end

%%
figure; 
for j = 1:length(unq_expgroups)
    exp_grp_j = unq_expgroups{j};
    inds = strcmp(pooled_table.ExpGroup, exp_grp_j);
    subplot(1,3,1); hold on; 
    plot(pooled_table.dV_TR_change(inds), pooled_table.per_nspk_change(inds), 'o', 'DisplayName', exp_grp_j);
    xlabel('dV change'); ylabel('\% nspk change');
    
    
    subplot(1,3,2); hold on; 
    plot(pooled_table.V_T_change(inds), pooled_table.per_nspk_change(inds), 'o', 'DisplayName', exp_grp_j);
    xlabel('$V_T$ change'); ylabel('\% nspk change');
    
    subplot(1,3,3); hold on;
    plot(pooled_table.V_R_change(inds), pooled_table.per_nspk_change(inds), 'o', 'DisplayName', exp_grp_j);
    xlabel('$V_R$ change'); ylabel('\% nspk change');
end
legend('show');

%%

[sorted_dVchange, sorted_ind] = sort(pooled_table.dV_TR_change);
relative_dVchange = pooled_table.dV_TR_change ./ pooled_table.dV_TR_base;

figure; hold on;
stem(sorted_dVchange, 'filled', '-ko', 'linewidth', 2, 'markersize', 5);
yyaxis right;
stem(100*relative_dVchange(sorted_ind), 'filled', '-ro')

%%
X_group = changed_table.ExpGroup; 
group_var = 'chol'; 

for i = 1:length(X_group)
    if contains(X_group{i}, group_var)
        X_group{i} = group_var;
    end
end
% 
% anovan(changed_table.n_spk_change, {X_group, changed_table.V_R_change, changed_table.V_T_change}, ...
%     'continuous', [2,3], ...
%     'model',1,'varnames',{'ExpGroup','V_R','V_T'});


anovan(changed_table.n_spk_change, {changed_table.ExpGroup, changed_table.V_R_change, changed_table.V_T_change}, ...
    'continuous', [2,3], ...
    'model',2,'varnames',{'IsChol','V_R','V_T'});

%% 

% anovan(pooled_table.n_spk_change, ...
%     {pooled_table.exp_chol, pooled_table.exp_elec, ...
%         pooled_table.dV_TR_change}, ...
%     'continuous', [3], ...
%     'model',1,'varnames',{'ExpChol','ExpElec','dV'});

% fitlm(pooled_table, 'n_spk_change ~ 1 + exp_chol + exp_elec + dV_TR_change')

%%

% X_group = changed_table.ExpGroup; 
% chol_group_ind = contains(X_group, 'chol');
% 
% [h,p] = ttest2(changed_table.V_T_change(chol_group_ind), changed_table.V_T_change(~chol_group_ind));
% [h,p] = ttest(changed_table.V_T_change(strcmpi(X_group, 'cholelec')));

%% HotellingT2 test
g1 = 'elec';
g2 = 'chol'; 

X_src = [pooled_table.dV_TR_change, pooled_table.n_spk_change];

ind_g1 = strcmp(pooled_table.ExpGroup,g1);
ind_g2 = strcmp(pooled_table.ExpGroup,g2);

group_ind = [ones(sum(ind_g1), 1); 2*ones(sum(ind_g2), 1)];

X = [group_ind, [X_src(ind_g1,:); X_src(ind_g2,:)]];
HotellingT2(X, 0.05);

%% Chi2 test of indepdendent

% effect_var = 'dV_TR_change';
% effect_cutoffs = 0.005:0.005:2;

effect_var = 'dV_TR_percent_change';
effect_cutoffs = 0.05:0.05:5;

% cutoff_type = 0;

figure; hold on;

for cutoff_type = [-1,0,1]
    p_vals = arrayfun(@(effect_cutoff) chi2test_groupeffect(pooled_table, effect_var, effect_cutoff, cutoff_type), effect_cutoffs);
    
    switch cutoff_type
        case 0
            lgnd = '(neg, 0, pos)';
        case -1
            lgnd = '(neg, nonneg)';
        case 1
            lgnd = '(nonpos, pos)';
    end
    plot(effect_cutoffs, p_vals, 'displayname', lgnd);
end

yline(0.05, '--k');
xlabel('cutoff');
ylabel('pval');
legend('show');

title(effect_var, 'Interpreter', 'none')

%%
% effect_var = 'dV_TR_change'; 
% effect_cutoff = 0.6;
% var_latex = '$\Delta V_{TR}^{\Delta}$';

effect_var = 'dV_TR_percent_change'; 
effect_cutoff = 3;
var_latex = '$\Delta V_{TR}^{\%\Delta}$';

[pval, chi2aux] = chi2test_groupeffect(pooled_table, effect_var, effect_cutoff, 0);

sign_group_data = chi2aux.data; 

pie_data = sign_group_data.Variables;
pie_titles = sign_group_data.Properties.RowNames;
pie_groups = sign_group_data.Properties.VariableNames;

figure; 
tiledlayout(1,length(pie_titles));
for i = 1:length(pie_titles)
    ax = nexttile;
    data_ith = pie_data(i,:);
    n_sample_ith = sum(pie_data(i,:));
    pie(ax, data_ith, '$%.1f\\%%$');
    set(findobj(ax,'type','text'), 'FontSize', 20);
    set(findobj(ax,'type','patch'), 'LineWidth', 5, 'EdgeColor', 'w');
    ax.Colormap = [0.05,0.05,0.05;0.4,0.4,0.4;0.9,0.9,0.9];
    title(sprintf('%s (n=%d)', pie_titles{i}, n_sample_ith));
end

lgd = legend(pie_groups);
lgd.Layout.Tile = 'east';

title(lgd, sprintf('split by \n %s = %.1f', var_latex, effect_cutoff));

%%  
    
    
function res = compare_expgroup_anov1(tbl, x_var, y_var, group_var)
X = tbl.(x_var); 
Y = tbl.(y_var); 

if strcmpi(group_var,'none') == 0
    for i = 1:length(X)
        if contains(X{i}, group_var) 
            X{i} = group_var;
        end
    end
end

[p,~,stats] = anova1(Y, X, 'off');
[c,~,~,gnames] = multcompare(stats);
expgroup_interaction = [gnames(c(:,1)), gnames(c(:,2)), num2cell(c(:,6))];

res = struct;
res.p_val_anova = p; 
res.expgroup_interaction = expgroup_interaction;

end

function [pval,aux] = chi2test_groupeffect(tbl, effect_var, effect_cutoff, cutoff_type)


[EG, EG_cat] = grp2idx(tbl.ExpGroup);
X = tbl.(effect_var);
SG = sign(X);

effect_cutoff = abs(effect_cutoff);

switch cutoff_type
    case 0
        SG(abs(X) < effect_cutoff) = 0;
        sign_group = {'neg','0','pos'};
    case 1
        SG(X < effect_cutoff) = -1;
        sign_group = {'nonpos','pos'};
    case -1
        SG(X > -effect_cutoff) = 1;
        sign_group = {'neg','nonneg'};
    otherwise
        error('Only -1,0,1 are accepted');     
end
   
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