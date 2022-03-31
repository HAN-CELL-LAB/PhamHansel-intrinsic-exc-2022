clc; clear; close all;
run start_up.m

graphic_setdefault(15, ...
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


%% Gor through and select post time 
% next_key = 'N'; 
% 
% tmp_xlines = gobjects(1,2); 
% 
% figure; 
% for i = 1:select_num_cells
%     sel_obj = select_analysis(i); 
%     sel_id = select_cellids(i); 
%     
%     subplot(212); hold on; 
%     plot(sel_obj.time_vec, sel_obj.Vthres_first, 'or', 'displayname', 'Vt');
%     plot(sel_obj.time_vec, sel_obj.Vrest_1, 'ok', 'displayname', 'Vr');
%     xline(0); 
%     legend('show');
%     subplot(211); 
%     hold on;
%     plot(sel_obj.time_vec, sel_obj.num_spikes, 'ok');
%     ylabel('n spk'); 
%     yyaxis right; 
%     plot(sel_obj.time_vec, sel_obj.Vthres_first - sel_obj.Vrest_1, 'or');
%     ylabel('Vt - Vr'); 
%     xline(0);
%     
%     title(sel_id); 
%     
%     cur_key = ''; 
%     
%     while ~strcmpi(cur_key, next_key)
%         delete(tmp_xlines); 
%         [t_selpost, ~] = ginput(2);         
%         t_selpost = min(max(round(t_selpost),1),sel_obj.time_vec(end))';
%         tmp_xlines = arrayfun(@(x) xline(x,'--k'), t_selpost, 'uni', 1);
%         title(sprintf('post = %d to %d', t_selpost(1), t_selpost(2)));
%         waitforbuttonpress
%         cur_key = get(gcf, 'currentcharacter');
%     end
%     
%     cell_selection.post{strcmp(cell_selection.cell_id, sel_id)} = t_selpost;
%     clf; 
% end
% 
% cell_selection.base = cellfun(@(x) sprintf('[%d,%d]',min(x),max(x)), cell_selection.base, 'uni', 0);
% cell_selection.post = cellfun(@(x) sprintf('[%d,%d]',min(x),max(x)), cell_selection.post, 'uni', 0);
% file_name = 'data/gill-data/Layer II.III S1 BC IP/processed/cell-selection.csv';
% writetable(cell_selection, file_name);

%% Select measures
% select_measures = {'num_spikes', 'Vthres_first', 'Vrest_1', 'Rin', 'fAHP_Vm_mean', 'dVthres_rest', 'dVfAHP_rest'}; 

select_measures = {'num_spikes', 'Vthres_first', 'Vrest_1', 'dVthres_rest', 'Rin'}; 

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
        tmp_struct.([measure_j '_change']) = post_j - base_j;
    end
    pooled_analysis{i} = tmp_struct;
end

pooled_analysis = structarray_to_struct(vertcat(pooled_analysis{:})); 
fields_pooled = fieldnames(pooled_analysis);
values_pooled = struct2cell(pooled_analysis);
%%
latex_changed = struct(...
    'num_spikes', 'n_{spk}', ...
    'Vthres_first', 'V_t', ...
    'Vrest_1', 'V_r', ...
    'Rin', 'R_i', ...
    'fAHP_Vm_mean', 'AHP', ...
    'dVthres_rest', '\Delta V_{t}', ...
    'dVfAHP_rest', '\Delta V_{AHP}');
latex_pooled = fields_pooled; 
latex_pooled_struct = structfun(@(x) '', pooled_analysis, 'uni', 0);
for i = 1:length(latex_pooled)
    raw_name = latex_pooled{i}; 
    main_name = strsplit(raw_name,{'_base','_change'});
    main_name = latex_changed.(main_name{1}); 
    if contains(raw_name, '_base') 
        latex_pooled{i} = sprintf('$%s^0$',main_name);
    elseif contains(raw_name, '_change') 
        latex_pooled{i} = sprintf('$%s^{\\mathbf{\\Delta}}$',main_name);
    end
    latex_pooled_struct.(raw_name) = latex_pooled{i};
end
% %% 
% [pooled_rhos, pooled_pval] = corrcoef(horzcat(values_pooled{:}));
% 
% figure; 
% subplot(121); hold on; 
% image_with_strict_limits(pooled_rhos);
% caxis([-1,1]); colormap(gca, return_colorbrewer('RdBu_r',100));
% daspect([1,1,1]);
% set(gca, ...
%     'xtick', 1:length(fields_pooled), 'xticklabel',latex_pooled,  ...
%     'XTickLabelRotation', 90, ...
%     'ytick', 1:length(fields_pooled), 'yticklabel',latex_pooled,  ...
%     'TickLabelInterpreter', 'latex')
% colorbar; 
% title('$\rho$ (pearson corr coeff)')
% subplot(122); hold on; 
% image_with_strict_limits(log10(pooled_pval));
% caxis([-2,0]); colormap(gca, return_colorbrewer('Blues_r',100));
% daspect([1,1,1]); colorbar; 
% set(gca, ...
%     'xtick', 1:length(fields_pooled), 'xticklabel',latex_pooled,  ...
%     'XTickLabelRotation', 90, ...
%     'ytick', 1:length(fields_pooled), 'yticklabel',latex_pooled,  ...
%     'TickLabelInterpreter', 'latex')
% title('$log_{10}$ (p-val)');
% export_fig('raw-rho-pval', '-dpng', '-r200', '-p0.02')

%%
% filt_rhos = triu(pooled_rhos .* (pooled_pval < 0.01));
% figure; hold on; 
% image_with_strict_limits(filt_rhos);
% caxis([-1,1]); colormap(gca, return_colorbrewer('RdBu_r',100));
% daspect([1,1,1]);
% set(gca, ...
%     'xtick', 1:length(fields_pooled), 'xticklabel',latex_pooled,  ...
%     'XTickLabelRotation', 90, ...
%     'ytick', 1:length(fields_pooled), 'yticklabel',latex_pooled,  ...
%     'TickLabelInterpreter', 'latex')
% colorbar; 
% title('$\rho$ (p-val $< 0.01$)')

% export_fig('rho-filt', '-dpng', '-r200', '-p0.02');
%%
% [rf, cf] = find(filt_rhos); 
% figure('defaultaxesfontsize', 15, ...
%         'DefaultAxesTitleFontSize', 1, ...
%     'DefaultAxesLabelFontSize', 1)
% for i = 1:length(rf)
%     subplot(3,5,i); hold on;
%     rn = fields_pooled{rf(i)};
%     cn = fields_pooled{cf(i)};
%     scatter(pooled_analysis.(rn), pooled_analysis.(cn), ...
%         200, [0.1,0.1,0.1], 'filled', 'markeredgecolor', 'none', ...
%     'markerfacealpha',0.3); 
%     xlabel(latex_pooled{rf(i)});
%     ylabel(latex_pooled{cf(i)});
%     axis square; 
% end
% 
% despline('all');

% export_fig('pair-filt', '-dpng', '-r200', '-p0.02');
%% 
nrows = 3;
ncols = 4; 

pval_signif = 0.05;
rest_splt = 1:(nrows * ncols);
rest_splt = rest_splt(mod(rest_splt,ncols) ~= 1);

for i = 1%:length(fields_pooled)
    figure; 
    
    field_name = fields_pooled{i};
    field_vp = regexprep(field_name, {'_base','_change'}, ''); 
    
    subplot(nrows,ncols,[1,5,9]); hold on; 
    
    yv = [pooled_analysis.([field_vp, '_base']), pooled_analysis.([field_vp, '_base'])+pooled_analysis.([field_vp, '_change'])];
    xv = repmat([1,2], [size(yv,1),1]);
    xv = xv + 0.2*(rand(size(xv)) - 0.5);
    vp = violinplot(yv, {'base (mean)','post (median)'}, ...
        'ShowData', false, 'ViolinColor', [0.7,0.7,0.7], ...
        'ShowNotches', true, 'BoxColor', [0.8,0.8,0.8], ...
        'EdgeColor', [1,1,1]);
    
    plot(xv',yv', 'color', [0.1,0.1,0.1,0.1], 'linewidth', 2);
    scatter(xv(:), yv(:), 50, [0.1,0.1,0.1], 'filled', 'markeredgecolor', 'none', ...
        'markerfacealpha',0.3);
    set(gca, 'xtick', [1,2], 'xticklabel', {'base', 'post'});
    ylabel(sprintf('$%s$', latex_changed.(field_vp)));
    despline;
    
    remain_fields = fields_pooled(~strcmp(fields_pooled,field_name)); 
    xv = pooled_analysis.(field_name); 
    xlbl = latex_pooled_struct.(field_name); 
    cnt_splt = 1; 
    for j = 1:length(remain_fields)
        field_j = remain_fields{j}; 
        subplot(nrows,ncols,rest_splt(cnt_splt)); hold on; 
        yv = pooled_analysis.(field_j);
        ylbl = latex_pooled_struct.(field_j);
        
        [cc,cpval] = corrcoef(xv,yv); cc = cc(1,2); cpval = cpval(1,2);
        cc_str = sprintf('$\\rho = %.3f$ and p-val $= %.3f$', cc, cpval);
        
        if cpval <= pval_signif
            highlight_color = [0.9,0,0];
        else
            highlight_color = [0.3,0.3,0.3];
        end
        
        line_fit = polyfit(xv, yv, 1);
        xFit = linspace(min(xv)-1, max(xv)+1, 100);
        yFit = polyval(line_fit , xFit);
        hold on;
        plot(xFit, yFit, '-', 'color', [0.2,0.2,0.2,0.1], 'LineWidth', 5);
        
        scatter(xv, yv, 100, [0.8,0.8,0.8], 'filled', 'markeredgecolor', 'none', ...
            'markerfacealpha', 0.9);
        ylabel(ylbl);
        if rest_splt(cnt_splt) > nrows*(ncols-1)
            xlabel(xlbl);
        end
        xline(0, '--', 'linewidth', 2);
        yline(0, '--', 'linewidth', 2);
        
        xlim([floor(min(xv)), ceil(max(xv))] + 0.1*[-1,1]*(max(xv)-min(xv)));
        ylim([floor(min(yv)), ceil(max(yv))] + 0.1*[-1,1]*(max(yv)-min(yv)));
        
       
        set(gca, 'xcolor', highlight_color, 'ycolor', highlight_color);
        title(cc_str, 'Color', highlight_color);
        
        cnt_splt = cnt_splt + 1; 
    end
    
%     export_fig(fullfile(fig_path, field_name),  '-dpng', '-r200', '-p0.02');
%     close; 
end

%%
% select_measures = {'num_spikes', 'latency', 'fAHP_Vm_mean', 'Vthres_first', 'Vrest_1', 'Rin'}; 
% 
% basetime_fun = @(x,y) mean(y(x >= -3 & x <= -1), 'omitnan');
% posttime_fun = @(x,y) median(y(x >= sum(x>0)*0.8 & x <= sum(x>0)*1.0), 'omitnan');
% 
% pooled_analysis = cellfun(@(s) ...
%     struct(s, arrayfun(@(x) struct(...
%     'base', basetime_fun(x.time_vec, x.(s)), ...
%     'post', posttime_fun(x.time_vec, x.(s)), ...
%     'change', posttime_fun(x.time_vec, x.(s)) - basetime_fun(x.time_vec, x.(s))), ...
%     select_analysis, 'uni', 1)), ...
%     select_measures, 'uni', 0);
% pooled_analysis = mergefield_struct(pooled_analysis{:});
% pooled_analysis = structfun(@structarray_to_struct, pooled_analysis, 'uni', 0);
% 
% % filt_sel = @(x) ...
% %     x.Vthres_first.base <= -35 & ...
% %     x.Vthres_first.post <= -35 & ...
% %     x.num_spikes.change > 0;
% filt_sel = @(x) x.Vthres_first.base  >= -1000; 
% sel_pool = filt_sel(pooled_analysis);
% pooled_analysis = structfun(@(x) structfun(@(y) y(sel_pool), x, 'uni', 0), pooled_analysis, 'uni', 0); 
% 
% 
% xv = pooled_analysis.Vthres_first.change; 
% yv = pooled_analysis.num_spikes.change; 
% [cc,cpval] = corrcoef(xv,yv); cc = cc(1,2); cpval = cpval(1,2);
% cc_str = sprintf('$\\rho = %.3f$\np-val $= %.3f$', cc, cpval);  
% cv = pooled_analysis.Rin.change; 
% rng_c = [-1,1]*max(abs(cv));
% cmap = return_colorbrewer('RdBu_r', 50)*0.9; 
% cv_mat = discretize_colorlevels(cv, rng_c, cmap); 
% 
% figure;
% subplot(1,3,[2,3]); hold on; 
% 
% line_fit = polyfit(xv, yv, 1);
% xFit = linspace(min(xv)-1, max(xv)+1, 100);
% yFit = polyval(line_fit , xFit);
% hold on;
% plot(xFit, yFit, '-', 'color', [0.2,0.2,0.2,0.1], 'LineWidth', 5);
% 
% scatter(xv, yv, 500, cv_mat, 'filled', 'markeredgecolor', 'none', ...
%     'markerfacealpha', 0.9); 
% xlabel('$ \Delta V_{\mathrm{thres}}^{(1)} (mV)$');
% ylabel('$\Delta$ \# spike'); 
% xline(0, '--', 'linewidth', 2);
% yline(0, '--', 'linewidth', 2);
% 
% xlim([floor(min(xv)), ceil(max(xv))]);
% ylim([floor(min(yv)), ceil(max(yv))]);
% 
% colormap(cmap); 
% caxis(rng_c); 
% cbar = colorbar; 
% xlabel(cbar, '$\Delta R_{\mathrm{in}}$', 'Interpreter', 'latex');
% cbar.Position = cbar.Position .* [1,1,0.5,0.5] + [0.08,0,0,0];
% 
% text(1,8, cc_str, 'FontSize', 30, 'Interpreter', 'latex')
% title('$\Delta = $ \textbf{post}(median, from 80\%) $-$ \textbf{base}(mean, last 3 traces)');
% 
% 
% axis square; 
% despline([2,5]); 
% 
% yv = [pooled_analysis.Vthres_first.base, pooled_analysis.Vthres_first.post]; 
% xv = repmat([1,2], [size(yv,1),1]);
% xv = xv + 0.1*(rand(size(xv)) - 0.5);
% 
% subplot(1,3,1); hold on;
% 
% vp = violinplot(yv, {'base (mean)','post (median)'}, ...
%     'ShowData', false, 'ViolinColor', [0.7,0.7,0.7], ...
%     'ShowNotches', true, 'BoxColor', [0.8,0.8,0.8], ...
%     'EdgeColor', [1,1,1]); 
% 
% plot(xv',yv', 'color', [0.1,0.1,0.1,0.1], 'linewidth', 5); 
% scatter(xv(:), yv(:), 200, [0.1,0.1,0.1], 'filled', 'markeredgecolor', 'none', ...
%     'markerfacealpha',0.3); 
% set(gca, 'xtick', [1,2], 'xticklabel', {'base', 'post'}); 
% ylabel('$V_{\mathrm{thres}}^{(1)} (mV)$');
% despline; 
% 
% export_fig('thres-1-vs-nspk-vs-Rin', '-dpng', '-r200', '-p0.02');
% close;
%%
