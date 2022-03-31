run start_up.m;
%% Load data
load('data/discrim_vs_recog_data.mat', ...
    'results_all', 'results_dim_description', ...
    'def_opts', 'param_opts');

fig_path = fullfile('figures/official');
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

%% Plot lines
common_norm_dthres_range = [-0.5, 0.5]; 

plot_fieldpairs = {'recognition_FPR', 'recognition_TPR', 'recognition_postradeoffs'};
title_name = struct; 
title_name.recognition_FPR = 'recog. FPR';
title_name.discrimination_accuracy = 'discrim. acc';
title_name.recognition_TPR = 'recog. TPR';
title_name.recognition_postradeoffs = 'TPR - FPR';

cbar_name = struct; 
cbar_name.discrimination_accuracy = 'discrim acc';
cbar_name.recognition_FPR = 'recog. FPR';
cbar_name.recognition_TPR = 'recog. TPR';
cbar_name.recognition_postradeoffs = 'TPR - FPR';

linestyles = {'-', '-', '-'}; 

cmap_pairs = {return_colorbrewer('Reds', length(percent_complete_input_vec)) * 0.9, ...
    return_colorbrewer('Blues', length(percent_complete_input_vec)) * 0.9, ...
    return_colorbrewer('Greens', length(percent_complete_input_vec)) * 0.95};

nearest_select_alpha = [0.25, 0.50, 0.6, 0.75, 1];
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

fig_description = sprintf('\\color[rgb]{%f, %f, %f}%s\\color{black} vs \\color[rgb]{%f, %f, %f}%s vs (\\color[rgb]{%f, %f, %f}%s)', ...
    cmap_pairs{1}(end,:), title_name.(plot_fieldpairs{1}), ...
    cmap_pairs{2}(end,:), title_name.(plot_fieldpairs{2}), ...
    cmap_pairs{3}(end,:), title_name.(plot_fieldpairs{3}));

annotation('textbox', 'String', fig_description, ...
    'FontSize', 25, 'Interpreter', 'tex', 'LineStyle', 'none', ...
    'Position', [0, 0.92, 0.8, 0.08], ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');

for j = 1:length(select_novrl_ind)
    for i = 1:length(select_alpha_ind)
        
        ax = subplot(fig_nrows,fig_ncols,cnt_splt); hold on;
        
        plot([0,0],[0,1], ':', 'linewidth', 2, 'color', 0.6*ones(1,3));
       
        for k = 1:length(percent_complete_input_vec)
            
            cellfun(@(x,c, lsty) ...
                plot(dthres_vec, squeeze(x(:,select_alpha_ind(i),k,select_novrl_ind(j))), ...
                'linewidth', 2.5, 'color', [c(k,:),0.9], 'linestyle', lsty), ...
                plot_resultlist, cmap_pairs, linestyles);
        end
        xlabel('$\Delta\theta$');
        
        if i == 1
            ylabel(sprintf(['\\textit{\\# overlap = %d}' newline 'task performance'], select_novrl_vec(j)));
        end
        
        if j == length(select_novrl_ind)
            xlabel({'$\Delta\theta$', sprintf('$\\alpha_W = %.2f$', select_alpha_vec(i))});
        else
            set(gca, 'xcolor', 'none');
        end
        
        cnt_splt = cnt_splt + 1;

    end
    
end

last_ax_pos = get(gca, 'Position');

linkaxes(findobj(gcf, 'type', 'axes'), 'xy');
for ax = findobj(gcf, 'type', 'axes')'
    ax.Position(1) = ax.Position(1) - 0.05;
end

despline('all');
        
xlim(common_norm_dthres_range);
ylim([-0.6,1]);

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


pseudo_ax_1 = axes('Position', last_ax_pos, 'Visible', 'off');
pseudo_ax_2 = axes('Position', last_ax_pos + [0,0.1,0,0], 'Visible', 'off');
pseudo_ax_3 = axes('Position', last_ax_pos + [0,0.25,0,0], 'Visible', 'off');

colormap(pseudo_ax_1, cmap_pairs{1});
cbar_1 = colorbar(pseudo_ax_1, 'north');
cbar_1.Position = cbar_1.Position .* [1,1,0.85,1] + [0,-0.05,0,0];
cbar_1.Position(1) = sum(last_ax_pos([1,3])) - 0.02;
title(cbar_1,sprintf('\\color[rgb]{%f, %f, %f} %s ', cmap_pairs{1}(end,:), cbar_name.(plot_fieldpairs{1})));
xlabel(cbar_1,'% input complete')
set(cbar_1, cbar_props{:});

colormap(pseudo_ax_2, cmap_pairs{2});
cbar_2 = colorbar(pseudo_ax_2, 'north');
cbar_2.Position([1,3,4]) = cbar_1.Position([1,3,4]);
title(cbar_2,sprintf('\\color[rgb]{%f, %f, %f} %s ', cmap_pairs{2}(end,:), cbar_name.(plot_fieldpairs{2})));
% xlabel(cbar_2,'% input complete')
set(cbar_2,  cbar_props{:});

colormap(pseudo_ax_3, cmap_pairs{3});
cbar_3 = colorbar(pseudo_ax_3, 'north');
cbar_3.Position([1,3,4]) = cbar_1.Position([1,3,4]);
title(cbar_3,sprintf('\\color[rgb]{%f, %f, %f} %s ', cmap_pairs{3}(end,:), cbar_name.(plot_fieldpairs{3})));
% xlabel(cbar_3,'% input complete')
set(cbar_3,  cbar_props{:});

fig_name = fullfile(fig_path, 'Fig5a');
export_fig(fig_name,  '-pdf', '-p0.02');
pause(0.5); close;

%% Plot delta_theta_optim, TRP, TPR-FPR, discacc
title_name = struct; 
title_name.recognition_accuracy = 'recognition TPR';
title_name.recognition_postradeoffs = 'recognition TPR - FPR';
title_name.discrimination_accuracy = 'discrimination accuracy';

flip_dtheta = true; 
diff_from_optim_acc = 0.001;
plot_fieldpairs = {'recognition_accuracy', 'recognition_postradeoffs', 'discrimination_accuracy'};

plot_resultlist = {results_all.(plot_fieldpairs{1}), ...
    results_all.(plot_fieldpairs{2}),  ...
    results_all.(plot_fieldpairs{3})};

nearest_select_alpha = 0.3:0.1:1;
nearest_select_novrl = [0,20];

cmap = flipud(bone(length(nearest_select_alpha))) * 0.85;

select_alpha_ind = arrayfun(@(x) find_nearest(alpha_W_vec, x, 'ind'), nearest_select_alpha);
select_novrl_ind = arrayfun(@(x) find_nearest(num_overlap_per_group_vec, x, 'ind'), nearest_select_novrl);

select_alpha_vec = alpha_W_vec(select_alpha_ind);
select_novrl_vec = num_overlap_per_group_vec(select_novrl_ind);

figure('units', 'normalized', 'Position', [0.05,0.05,1,0.9], 'defaultaxesfontsize', 20);

cnt_splt = 1;
for i = 1:length(select_novrl_ind)
    for k = 1:length(plot_fieldpairs)
        subplot(length(select_novrl_ind),length(plot_fieldpairs),cnt_splt); hold on;
        
        plot([0,1],[0,0], '--', 'linewidth', 2, 'color', 0.4*ones(1,3));
        
        for j = 1:length(select_alpha_ind)
            acc_select = squeeze(plot_resultlist{k}(:,select_alpha_ind(j),:,select_novrl_ind(i)));
            
            diff_from_max_acc = abs(acc_select -  max(acc_select,[],1));
            
            dthres_discrim_max = arrayfun(@(x) ...
                dthres_vec(diff_from_max_acc(:,x) < diff_from_optim_acc), ...
                1:size(diff_from_max_acc,2), 'uni', 0);
            dthres_optim = cellfun(@(x) find_nearest(x, 0), dthres_discrim_max);
            
            plot(percent_complete_input_vec, dthres_optim, 'linewidth', 3, ...
                'color', [cmap(j,:), 0.7]);
        end
        
        
        xlabel('\% complete');
        if k == 1
            ylabel({sprintf('\\# overlap = %d', select_novrl_vec(i)), '$\Delta\theta_{\mathrm{optim}}$'});
        end
        
        if flip_dtheta
            ylabel('$\Delta\theta_{\mathrm{optim}}$');
            if k == 1                
                xlabel({sprintf('\\textit{\\# overlap = %d}', select_novrl_vec(i)), '\% complete'});
            end
        end
        
        if i == 1
            title(title_name.(plot_fieldpairs{k}));
        end    
        if flip_dtheta
            view([90 -90]);
        end
        
        cnt_splt = cnt_splt + 1;
        
    end
    
end

linkaxes(findall(gcf, 'type', 'axes'), 'xy');
ylim([-0.8,0.2]); 
xlim([0,1]);

despline('all',[0.5,0.3]);
    
subplot(length(select_novrl_ind),length(plot_fieldpairs),6);
ax_pos = get(gca, 'Position');
cbar = colorbar(gca);

select_cbartickpos = 1:2:length(select_alpha_vec);
cbar_tickval = arrayfun(@(x) sprintf('%.1g', x), select_alpha_vec, 'uni', 0);
cbar_numtick = length(cbar_tickval);
cbar_tickpos = linspace(0.5/cbar_numtick,1-0.5/cbar_numtick,cbar_numtick);
cbar_props = {...
    'linewidth', 1.5, ...
    'ticklength', 0.015, ...
    'fontsize', 28, ...
    'ticks', cbar_tickpos(select_cbartickpos), ...
    'ticklabels', cbar_tickval(select_cbartickpos)};

colormap(gca, cmap);
cbar.Position(1) = ax_pos(1) + ax_pos(3) - 0.04;
cbar.Position = cbar.Position .* [1,1,1.25,0.65] + [0.06,0.04,0,0];
title(cbar, '$\alpha_W$', 'interpreter', 'latex', 'fontsize', 35);
set(cbar, cbar_props{:});

fig_name = fullfile(fig_path, 'Fig6a');
export_fig(fig_name,  '-pdf', '-p0.02');
pause(0.5); close;
