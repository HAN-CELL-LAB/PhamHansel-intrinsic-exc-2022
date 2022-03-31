run start_up.m;
%% Load data
load('data/discrim_vs_recog_data.mat', ...
    'results_all', 'results_dim_description', ...
    'def_opts', 'param_opts');

fig_path = fullfile('figures/discrim-vs-recog/inprogress');
if ~exist(fig_path, 'dir')
    mkdir(fig_path);
end

%% Adding L2 norm of recognition rates
results_all.recognition_L2_truerates = sqrt(results_all.recognition_TPR.^2 + results_all.recognition_TNR.^2);
results_all.recognition_L2_falserates = sqrt(results_all.recognition_FPR.^2 + results_all.recognition_FNR.^2);

results_all.recognition_L1_truerates = (results_all.recognition_TPR + results_all.recognition_TNR);
results_all.recognition_L1_falserates = (results_all.recognition_FPR + results_all.recognition_FNR);

tradeoff_lambda = 1; 
results_all.recognition_postradeoffs = results_all.recognition_TPR - tradeoff_lambda * results_all.recognition_FPR; 

%% Vectors of parameters
dthres_vec = linspace(def_opts.rng_dthres(1), def_opts.rng_dthres(2), def_opts.num_dthres);
alpha_W_vec = param_opts.alpha_W;
percent_complete_input_vec = param_opts.percent_complete_input;
num_overlap_per_group_vec = param_opts.num_overlap_per_group;

plot_fields = fieldnames(results_all);

%% Plot ROC lines
common_norm_dthres_range = [-0.5, 0.5]; 

plot_fieldpairs = {'recognition_FPR', 'recognition_TPR'};
title_name = struct; 
title_name.recognition_FPR = 'recog. FPR';
title_name.recognition_TPR = 'recog. TPR';

cbar_name = struct; 
cbar_name.recognition_FPR = 'recog. FPR';
cbar_name.recognition_TPR = 'recog. TPR';

linestyles = {'-', '-'}; 

cmap_pairs = {return_colorbrewer('Reds', length(percent_complete_input_vec)) * 0.9, ...
    return_colorbrewer('Blues', length(percent_complete_input_vec)) * 0.9};

nearest_select_alpha = [0.25, 0.50, 0.6, 0.75, 1];
% nearest_select_alpha = [0.50, 0.6, 0.75, 0.85, 1];
nearest_select_novrl = [0, 20];

select_alpha_ind = arrayfun(@(x) find_nearest(alpha_W_vec, x, 'ind'), nearest_select_alpha);
select_novrl_ind = arrayfun(@(x) find_nearest(num_overlap_per_group_vec, x, 'ind'), nearest_select_novrl);

select_alpha_vec = alpha_W_vec(select_alpha_ind);
select_novrl_vec = num_overlap_per_group_vec(select_novrl_ind);

fig_ncols = length(select_alpha_ind);
fig_nrows = length(select_novrl_vec);

plot_resultlist = cellfun(@(x) results_all.(x), plot_fieldpairs, 'uni', 0);

figure('units', 'normalized', ...
    'Position', [0.05,0.05,1,0.65], ...
    'DefaultAxesFontSize', 18);
cnt_splt = 1;

fig_description = 'ROC for recognition';
annotation('textbox', 'String', fig_description, ...
    'FontSize', 25, 'Interpreter', 'tex', 'LineStyle', 'none', ...
    'Position', [0, 0.9, 1, 0.08], ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');


cmap = flipud(bone(length(percent_complete_input_vec))) * 0.85;
for j = 1:length(select_novrl_ind)
    for i = 1:length(select_alpha_ind)
        
        ax = subplot(fig_nrows,fig_ncols,cnt_splt); hold on;
        
        plot([0,1],[0,1], '--', 'linewidth', 2, 'color', 'k');
       
        for k = 1:length(percent_complete_input_vec)
            x = plot_resultlist{1}(:,select_alpha_ind(i),k,select_novrl_ind(j));
            y = plot_resultlist{2}(:,select_alpha_ind(i),k,select_novrl_ind(j));
            p_cmp = cmap(k,:);
            plot(x, y, 'color', p_cmp, 'linewidth', 3);
        end
        xlabel('FPR');
        
        ax.Position(2) = ax.Position(2) - 0.1*(1.5-j); 
        
        if i == 1
            ylabel(sprintf(['\\textit{\\# overlap = %d}' newline 'TPR'], select_novrl_vec(j)));
        end
        
        if j == length(select_novrl_ind)
            xlabel({'FPR', sprintf('$\\alpha_W = %.2f$', select_alpha_vec(i))});
        else
            set(gca, 'xcolor', 'none');
        end
        
        cnt_splt = cnt_splt + 1;

    end
    
end

despline('all');
linkaxes(findall(gcf, 'type', 'axes'), 'xy')


last_ax_pos = get(gca, 'Position');

pseudo_ax_1 = axes('Position', last_ax_pos, 'Visible', 'off');

select_cbartickpos = [1,4,7,10];

cbar_tickval = percent_complete_input_vec;
cbar_numtick = length(cbar_tickval);
cbar_tickpos = linspace(0.5/cbar_numtick,1-0.5/cbar_numtick,cbar_numtick);
cbar_props = {...
    'linewidth', 1.5, ...
    'ticklength', 0.015, ...
    'fontsize', 18, ...
    'ticks', cbar_tickpos(select_cbartickpos), ...
    'ticklabels', cbar_tickval(select_cbartickpos)};


colormap(pseudo_ax_1, cmap);
cbar_1 = colorbar(pseudo_ax_1);
cbar_1.Position = cbar_1.Position .* [1,1,0.9,0.75] + [0,0.05,0,0];
cbar_1.Position(1) = sum(last_ax_pos([1,3])) + 0.01;
xlabel(cbar_1,'% input complete')
set(cbar_1, cbar_props{:});

fig_name = fullfile(fig_path, 't4');
export_fig(fig_name,  '-r300', '-p0.02');
pause(0.5); close;


%%
common_norm_dthres_range = [-0.5, 0.5]; 

plot_fieldpairs = {'recognition_FPR', 'recognition_TPR'};
title_name = struct; 
title_name.recognition_FPR = 'recog. FPR';
title_name.recognition_TPR = 'recog. TPR';

cbar_name = struct; 
cbar_name.recognition_FPR = 'recog. FPR';
cbar_name.recognition_TPR = 'recog. TPR';

linestyles = {'-', '-'}; 

cmap_pairs = {return_colorbrewer('Reds', length(percent_complete_input_vec)) * 0.9, ...
    return_colorbrewer('Blues', length(percent_complete_input_vec)) * 0.9};

nearest_select_alpha = [0.25, 0.50, 0.6, 0.75, 1];
% nearest_select_alpha = [0.50, 0.6, 0.75, 0.85, 1];
nearest_select_novrl = [0, 20];

select_alpha_ind = arrayfun(@(x) find_nearest(alpha_W_vec, x, 'ind'), nearest_select_alpha);
select_novrl_ind = arrayfun(@(x) find_nearest(num_overlap_per_group_vec, x, 'ind'), nearest_select_novrl);

select_alpha_vec = alpha_W_vec(select_alpha_ind);
select_novrl_vec = num_overlap_per_group_vec(select_novrl_ind);

fig_ncols = length(select_alpha_ind);
fig_nrows = length(select_novrl_vec);

plot_resultlist = cellfun(@(x) results_all.(x), plot_fieldpairs, 'uni', 0);

figure('units', 'normalized', ...
    'Position', [0.05,0.05,1,0.65], ...
    'DefaultAxesFontSize', 18);
cnt_splt = 1;

% fig_description = sprintf('\\color[rgb]{%f, %f, %f}%s\\color{black} vs \\color[rgb]{%f, %f, %f}%s)', ...
%     cmap_pairs{1}(end,:), title_name.(plot_fieldpairs{1}), ...
%     cmap_pairs{2}(end,:), title_name.(plot_fieldpairs{2}));
% 
% annotation('textbox', 'String', fig_description, ...
%     'FontSize', 25, 'Interpreter', 'tex', 'LineStyle', 'none', ...
%     'Position', [0, 0.92, 0.8, 0.08], ...
%     'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');

select_dthres_vec_ind = 12:(50-11);
select_dthres = dthres_vec(select_dthres_vec_ind);
cmap = return_colorbrewer('RdBu', length(select_dthres_vec_ind)) * 0.85;
for j = 1:length(select_novrl_ind)
    for i = 1:length(select_alpha_ind)
        
        ax = subplot(fig_nrows,fig_ncols,cnt_splt); hold on;
        
        plot([0,1],[0,1], '--', 'linewidth', 2, 'color', 'k');
       
        for k = 1:length(select_dthres_vec_ind)
            k_ind = select_dthres_vec_ind(k);
            x = squeeze(plot_resultlist{1}(k_ind,select_alpha_ind(i),:,select_novrl_ind(j)));
            y = squeeze(plot_resultlist{2}(k_ind,select_alpha_ind(i),:,select_novrl_ind(j)));
            p_cmp = cmap(k,:);
            plot(x, y, 'color', p_cmp, 'linewidth', 3);
        end
        xlabel('FPR');
        
        if i == 1
            ylabel(sprintf(['\\textit{\\# overlap = %d}' newline 'TPR'], select_novrl_vec(j)));
        end
        
        if j == length(select_novrl_ind)
            xlabel({'FPR', sprintf('$\\alpha_W = %.2f$', select_alpha_vec(i))});
        else
            set(gca, 'xcolor', 'none');
        end
        
        cnt_splt = cnt_splt + 1;

    end
    
end

despline('all');
linkaxes(findall(gcf, 'type', 'axes'), 'xy')


last_ax_pos = get(gca, 'Position');

pseudo_ax_1 = axes('Position', last_ax_pos, 'Visible', 'off');

select_cbartickpos = [1,4,7,10];

cbar_tickval = percent_complete_input_vec;
cbar_numtick = length(cbar_tickval);
cbar_tickpos = linspace(0.5/cbar_numtick,1-0.5/cbar_numtick,cbar_numtick);
cbar_props = {...
    'linewidth', 1.5, ...
    'ticklength', 0.015, ...
    'fontsize', 18};


colormap(pseudo_ax_1, cmap);
cbar_1 = colorbar(pseudo_ax_1);
cbar_1.Position = cbar_1.Position .* [1,1,0.9,0.75] + [0,0.05,0,0];
cbar_1.Position(1) = sum(last_ax_pos([1,3])) + 0.01;
xlabel(cbar_1,'$\Delta \theta$', 'Interpreter', 'latex')
caxis(select_dthres([1,end]))
set(cbar_1, cbar_props{:});

% fig_name = fullfile(fig_path, 't5');
% export_fig(fig_name,  '-r300', '-p0.02');
% pause(0.5); close;


%% Plot auc roc 
chance_alphaW = 1/3; 
chance_recogclassif = 1/2; 

plot_fieldpairs = {'recognition_FPR', 'recognition_TPR'};

plot_resultlist = {results_all.(plot_fieldpairs{1}), ...
    results_all.(plot_fieldpairs{2})};

cmap = flipud(bone(length(percent_complete_input_vec))) * 0.85;
select_dthres_ind = find(abs(dthres_vec) <= 0.52); 

figure('units', 'normalized', 'Position', [0.05,0.05,0.9,0.9], 'defaultaxesfontsize', 15);

cnt_splt = 1;
for ind_novrl = 1:(length(num_overlap_per_group_vec)-1)
    novrl = num_overlap_per_group_vec(ind_novrl); 
    
    subplot(3,3,ind_novrl); hold on; 
    for ind_inp = 1:length(percent_complete_input_vec)
        
        auc_vec = nan(1, length(alpha_W_vec));
        for ind_alpha = 1:length(alpha_W_vec)
            FPRv = squeeze(plot_resultlist{1}(select_dthres_ind,ind_alpha,ind_inp,ind_novrl));
            TPRv = squeeze(plot_resultlist{2}(select_dthres_ind,ind_alpha,ind_inp,ind_novrl));
            FPRv = sort(FPRv);
            TPRv = sort(TPRv);
            FPRv = [0;FPRv;1];
            TPRv = [0;TPRv;1];
            auc_vec(ind_alpha) = trapz(FPRv,TPRv);
        end
        plot(alpha_W_vec, auc_vec, 'linewidth', 3, 'color', cmap(ind_inp, :)); 
        
    end
    linechance_alphaW = xline(chance_alphaW, ':k', 'linewidth',2);
    linechance_recogclassif = yline(chance_recogclassif, '--k', 'linewidth',2);
    
    if ind_novrl == 1
        xlabel('$\alpha_W$');
        ylabel('recog. AUC ROC');
        legend([linechance_alphaW, linechance_recogclassif], ...
            '$\alpha_W=1/N_Y=1/3$', 'recog. chance');
        
    end
    
    title(sprintf('\\# overlap = %d',novrl));
end

despline('all');
linkaxes(findall(gcf,'type','axes'),'xy');


fig_name = fullfile(fig_path, 't6');
export_fig(fig_name,  '-r300', '-p0.02');
pause(0.5); close;