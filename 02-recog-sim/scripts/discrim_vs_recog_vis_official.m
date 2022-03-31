run start_up.m;
%% Load data
load('data/discrim_vs_recog_data.mat', ...
    'results_all', 'results_dim_description', ...
    'def_opts', 'param_opts');

fig_path = fullfile('figures/discrim-vs-recog');
if ~exist(fig_path, 'dir')
    mkdir(fig_path);
end

%% Vectors of parameters
dthres_vec = linspace(def_opts.rng_dthres(1), def_opts.rng_dthres(2), def_opts.num_dthres);
alpha_W_vec = param_opts.alpha_W;
percent_complete_input_vec = param_opts.percent_complete_input;
num_overlap_per_group_vec = param_opts.num_overlap_per_group;

plot_fields = fieldnames(results_all);

%% Plot lines

plot_fieldpairs = {'discrimination_accuracy', 'recognition_accuracy'};

cmap_pairs = {return_colorbrewer('Reds', length(percent_complete_input_vec)) * 0.9, ...
    return_colorbrewer('Blues', length(percent_complete_input_vec)) * 0.9};

nearest_select_alpha = [0.25, 0.50, 0.75, 1];
nearest_select_novrl = [0, 20];

select_alpha_ind = arrayfun(@(x) find_nearest(alpha_W_vec, x, 'ind'), nearest_select_alpha);
select_novrl_ind = arrayfun(@(x) find_nearest(num_overlap_per_group_vec, x, 'ind'), nearest_select_novrl);

select_alpha_vec = alpha_W_vec(select_alpha_ind);
select_novrl_vec = num_overlap_per_group_vec(select_novrl_ind);

fig_ncols = length(select_alpha_ind);
fig_nrows = length(select_novrl_vec);

plot_resultlist = {results_all.(plot_fieldpairs{1}), ...
    results_all.(plot_fieldpairs{2})};

figure('units', 'normalized', 'Position', [0.05,0.05,0.9,0.6]);
cnt_splt = 1;

fig_description = sprintf('\\color[rgb]{%f, %f, %f}%s\\color{black} vs \\color[rgb]{%f, %f, %f}%s', ...
    cmap_pairs{1}(end,:), regexprep(plot_fieldpairs{1}, '_', '-'), ...
    cmap_pairs{2}(end,:), regexprep(plot_fieldpairs{2}, '_', '-'));

annotation('textbox', 'String', fig_description, ...
    'FontSize', 20, 'Interpreter', 'tex', 'LineStyle', 'none', ...
    'Position', [0, 0.92, 1, 0.08], ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');

for j = 1:length(select_novrl_ind)
    for i = 1:length(select_alpha_ind)
        
        subplot(fig_nrows,fig_ncols,cnt_splt); hold on;
        
        plot([0,0],[0,1], '--', 'linewidth', 2, 'color', 0.4*ones(1,3));
        
        for k = 1:length(percent_complete_input_vec)
            
            cellfun(@(x,c) plot(dthres_vec, squeeze(x(:,select_alpha_ind(i),k,select_novrl_ind(j))), ...
                'linewidth', 3, 'color', c(k,:)), plot_resultlist, cmap_pairs);
        end
        xlabel('$\Delta\theta$');
        
        if i == 1
            ylabel(sprintf('\\# overlap = %d', select_novrl_vec(j)));
        end
        
        if j == length(select_novrl_ind)
            xlabel({'$\Delta\theta$', sprintf('$\\alpha_W = %.2f$', select_alpha_vec(i))});
        else
            set(gca, 'xcolor', 'none');
        end
        
        despline;
        
        cnt_splt = cnt_splt + 1;
    end
    
    
end

select_cbartickpos = [1,4,7,10];

cbar_tickval = percent_complete_input_vec;
cbar_numtick = length(cbar_tickval);
cbar_tickpos = linspace(0.5/cbar_numtick,1-0.5/cbar_numtick,cbar_numtick);
cbar_props = {...
    'fontsize', 15, ...
    'ticks', cbar_tickpos(select_cbartickpos), ...
    'ticklabels', cbar_tickval(select_cbartickpos)};

last_ax_pos = get(gca, 'Position');
pseudo_ax_1 = axes('Position', last_ax_pos, 'Visible', 'off');
pseudo_ax_2 = axes('Position', last_ax_pos, 'Visible', 'off');

colormap(pseudo_ax_1, cmap_pairs{1});
cbar_1 = colorbar(pseudo_ax_1);
cbar_1.Position = cbar_1.Position .* [1,1,1,0.8] + [0,0,0,0];
cbar_1.Position(1) = sum(last_ax_pos([1,3])) + 0.02;
ylabel(cbar_1, ['% complete' newline sprintf('(%s)', regexprep(plot_fieldpairs{1}, '_', '-'))]);
set(cbar_1, cbar_props{:});

colormap(pseudo_ax_2, cmap_pairs{2});
cbar_2 = colorbar(pseudo_ax_2);
cbar_2.Position = cbar_1.Position + [0,0.4,0,0];
ylabel(cbar_2, ['% complete' newline sprintf('(%s)', regexprep(plot_fieldpairs{2}, '_', '-'))]);
set(cbar_2,  cbar_props{:});

fig_name = fullfile(fig_path, sprintf('OFFICIAL-%s-AND-%s-lineplotOfalphaWAndnOverlap', plot_fieldpairs{1}, plot_fieldpairs{2}));
export_fig(fig_name,  '-r300', '-p0.02');
pause(0.5); close;


%% Plot delta_theta_optim
title_name = struct; 
title_name.discrimination_accuracy = 'discrimination accuracy';
title_name.recognition_accuracy = 'recognition TPR';

flip_dtheta = true; 
diff_from_optim_acc = 0.01;
plot_fieldpairs = {'discrimination_accuracy', 'recognition_accuracy'};

plot_resultlist = {results_all.(plot_fieldpairs{1}), ...
    results_all.(plot_fieldpairs{2})};

nearest_select_alpha = 0.2:0.1:1;
nearest_select_novrl = [0,20];

cmap = flipud(bone(length(nearest_select_alpha))) * 0.85;

select_alpha_ind = arrayfun(@(x) find_nearest(alpha_W_vec, x, 'ind'), nearest_select_alpha);
select_novrl_ind = arrayfun(@(x) find_nearest(num_overlap_per_group_vec, x, 'ind'), nearest_select_novrl);

select_alpha_vec = alpha_W_vec(select_alpha_ind);
select_novrl_vec = num_overlap_per_group_vec(select_novrl_ind);

figure('units', 'normalized', 'Position', [0.05,0.05,0.75,0.9], 'defaultaxesfontsize', 25);

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
ylim([-0.8,0.2]); xlim([0,1]);

despline('all',[0.5,0.3]);
    
subplot(length(select_novrl_ind),length(plot_fieldpairs),4);
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
cbar.Position = cbar.Position .* [1,1,1.2,0.7] + [0,0.05,0,0];
title(cbar, '$\alpha_W$', 'interpreter', 'latex');
set(cbar, cbar_props{:});

fig_name = fullfile(fig_path, sprintf('OFFICIAL-%s-AND-%s-lineplotOfalphaWAndnOverlap-dThetaOptim', plot_fieldpairs{1}, plot_fieldpairs{2}));
export_fig(fig_name,  '-r300', '-p0.02');
pause(0.5); close;

%%
diff_from_optim_acc = 0.01;
discrim_acc = results_all.discrimination_accuracy;

nearest_select_alpha = 0.2:0.1:1;
nearest_select_novrl = [0,20,40];

cmap = flipud(bone(length(nearest_select_alpha))) * 0.85;

select_alpha_ind = arrayfun(@(x) find_nearest(alpha_W_vec, x, 'ind'), nearest_select_alpha);
select_novrl_ind = arrayfun(@(x) find_nearest(num_overlap_per_group_vec, x, 'ind'), nearest_select_novrl);

select_alpha_vec = alpha_W_vec(select_alpha_ind);
select_novrl_vec = num_overlap_per_group_vec(select_novrl_ind);

figure('units', 'normalized', 'Position', [0.05,0.05,0.5,0.8], 'defaultaxesfontsize', 25);
hold on; linestyles = {'-', '--', ':', '-.'};

lgnd_objs = gobjects(length(select_novrl_ind),1);

cnt_splt = 1;
for i = 1:length(select_novrl_ind)
%     subplot(length(select_novrl_ind),1,cnt_splt); hold on;
    
    for j = 1:length(select_alpha_ind)
        acc_select = squeeze(discrim_acc(:,select_alpha_ind(j),:,select_novrl_ind(i)));
        max_select = max(acc_select,[],1);
        plt_handle = plot(percent_complete_input_vec, max_select, 'linewidth', 3, ...
            'color', cmap(j,:), 'linestyle', linestyles{i});
        if j == length(select_alpha_ind)
            lgnd_objs(i) = plt_handle;
            plt_handle.DisplayName = sprintf('\\# overlap = %d', select_novrl_vec(i));
        end
    end
    
    xlabel('\% complete');
    ylabel('accuracy');
    if i == 1
        title('max discrimination-accuracy');
    end
    
    cnt_splt = cnt_splt + 1;
    
end

xlim([0.1,1]);
despline('all',[1,10]);
% linkaxes(findall(gcf, 'type', 'axes'), 'xy');
% ax_pos = get(gca, 'Position');
% cbar = colorbar(gca);
% caxis(select_alpha_vec([1,end]));
% colormap(gca, cmap);
% cbar.Position(1) = ax_pos(1) + ax_pos(3) + 0.02;
% cbar.Position = cbar.Position .* [1,1,1,1];
% cbar.FontSize = 15;

ylim([-0.05,1]); 
legend(lgnd_objs, 'fontsize', 25);
% cbar = colorbar(gca);
% caxis(select_alpha_vec([1,end]));
% colormap(gca, cmap);
% cbar.Position(1) = ax_pos(1) + ax_pos(3) + 0.02;
% cbar.Position = cbar.Position .* [1,1,1,0.5];
% cbar.FontSize = 20;
% title(cbar, '$\alpha_W$', 'interpreter', 'latex');

fig_name = fullfile(fig_path, 'OFFICIAL-discrimination-max-lineplotOfalphaWAndnOverlap');
export_fig(fig_name,  '-r300', '-p0.02');
pause(0.5); close;
