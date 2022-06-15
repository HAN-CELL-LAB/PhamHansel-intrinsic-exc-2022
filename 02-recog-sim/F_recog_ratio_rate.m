clc; clear; close all;
run start_up.m;

%% Load data
load('data/recog_data.mat', ...
    'results_all', 'results_dim_description', ...
    'def_opts', 'param_opts');

fig_path = fullfile('figures/official');
if ~exist(fig_path, 'dir')
    mkdir(fig_path);
end


%% Vectors of parameters
dthres_vec = linspace(def_opts.rng_dthres(1), def_opts.rng_dthres(2), def_opts.num_dthres);
alpha_W_vec = param_opts.alpha_W;
percent_complete_input_vec = param_opts.percent_complete_input;
num_overlap_per_group_vec = param_opts.num_overlap_per_group;

plot_fields = fieldnames(results_all);

results_all.recognition_postradeoffs = results_all.recognition_TPR - results_all.recognition_FPR;
results_all.recognition_PLLR = log10((results_all.recognition_TPR + 1e-5) ./ (results_all.recognition_FPR + 1e-5));
results_all.recognition_nPRratio = results_all.recognition_TPR ./ (results_all.recognition_FPR + results_all.recognition_TPR);

%% Plot lines
common_norm_dthres_range = [-0.5, 0.5];

plot_fieldpairs_set = {...
    {'recognition_PLLR', 'recognition_TPR'},...
    {'recognition_nPRratio', 'recognition_TPR'}
    };

title_name = struct;
title_name.recognition_postradeoffs = 'recog. TPR - FPR';
title_name.recognition_TPR = 'recog. TPR';
title_name.recognition_PLLR = 'recog. log_{10}(TPR / FPR)';
title_name.recognition_nPRratio = 'recog. TPR / (TPR + FPR)';

cbar_name = struct;
cbar_name.recognition_postradeoffs = 'recog. TPR - FPR';
cbar_name.recognition_TPR = 'recog. TPR';
cbar_name.recognition_PLLR = 'recog. log_{10}(TPR / FPR)';
cbar_name.recognition_nPRratio = 'recog. TPR / (TPR + FPR)';

cmap_pairs = {return_colorbrewer('Reds', length(percent_complete_input_vec)) * 0.9, ...
    return_colorbrewer('Blues', length(percent_complete_input_vec)) * 0.9};

nearest_select_alpha = [0.25, 0.50, 0.6, 0.75, 1];
nearest_select_novrl = [0, 20];

select_alpha_ind = arrayfun(@(x) find_nearest(alpha_W_vec, x, 'ind'), nearest_select_alpha);
select_novrl_ind = arrayfun(@(x) find_nearest(num_overlap_per_group_vec, x, 'ind'), nearest_select_novrl);

select_alpha_vec = alpha_W_vec(select_alpha_ind);
select_novrl_vec = num_overlap_per_group_vec(select_novrl_ind);

fig_ncols = length(select_alpha_ind);
fig_nrows = length(select_novrl_vec);

for i = 1:length(plot_fieldpairs_set)
    plot_fieldpairs = plot_fieldpairs_set{i};
    plot_resultlist = {results_all.(plot_fieldpairs{1}), ...
        results_all.(plot_fieldpairs{2})};
    
    figure('units', 'normalized', ...
        'Position', [0,0.05,1,0.6], ...
        'DefaultAxesFontSize', 16);
    
    cnt_splt = 1;
    
    fig_description = sprintf('\\color[rgb]{%f, %f, %f}%s\\color{black} vs \\color[rgb]{%f, %f, %f}%s', ...
        cmap_pairs{1}(end,:), title_name.(plot_fieldpairs{1}), ...
        cmap_pairs{2}(end,:), title_name.(plot_fieldpairs{2}));
    
    annotation('textbox', 'String', fig_description, ...
        'FontSize', 25, 'Interpreter', 'tex', 'LineStyle', 'none', ...
        'Position', [0, 0.90, 1, 0.1], ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
    
    for j = 1:length(select_novrl_ind)
        for i = 1:length(select_alpha_ind)
            
            subplot(fig_nrows,fig_ncols,cnt_splt); hold on;
            
            plot([0,0],[0,1], ':', 'linewidth', 2, 'color', 0.6*ones(1,3));
            
            for k = 1:length(percent_complete_input_vec)
                cellfun(@(x,c) plot(dthres_vec, squeeze(x(:,select_alpha_ind(i),k,select_novrl_ind(j))), ...
                    'linewidth', 3, 'color', c(k,:)), plot_resultlist, cmap_pairs);
            end
            xlabel('$\Delta\theta$');
            
            if i == 1
                ylabel(sprintf(['\\textit{\\# overlap = %d}' newline 'task performance'], select_novrl_vec(j)));
            end
            
            if j == length(select_novrl_ind)
                xlabel('$\Delta\theta$');
            else
                set(gca, 'xcolor', 'none');
            end
            
            if j == length(select_novrl_ind)
                xlabel({'$\Delta\theta$', sprintf('$\\alpha_W = %.2f$', select_alpha_vec(i))});
            else
                set(gca, 'xcolor', 'none');
            end
            
            xlim(common_norm_dthres_range);
            despline;
            
            ax_pos = get(gca, 'Position');
            ax_pos = ax_pos + [-0.05,(j-1)*0.1-0.05,0,0];
            set(gca, 'Position', ax_pos);
            
            cnt_splt = cnt_splt + 1;
        end
        
        
    end
    
    linkaxes(findobj(gcf, 'type', 'axes'), 'xy');
    
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
    cbar_1.Position(1) = sum(last_ax_pos([1,3])) + 0.05;
    title(cbar_1,sprintf('\\color[rgb]{%f, %f, %f} %s ', cmap_pairs{1}(end,:), cbar_name.(plot_fieldpairs{1})));
    xlabel(cbar_1,'% input complete')
    set(cbar_1, cbar_props{:});
    
    colormap(pseudo_ax_2, cmap_pairs{2});
    cbar_2 = colorbar(pseudo_ax_2);
    cbar_2.Position = cbar_1.Position + [0,0.4,0,0];
    title(cbar_2,sprintf('\\color[rgb]{%f, %f, %f} %s ', cmap_pairs{2}(end,:), cbar_name.(plot_fieldpairs{2})));
    xlabel(cbar_2,'% input complete')
    set(cbar_2,  cbar_props{:});
    
    fig_name = sprintf('supp_%s-vs-%s.pdf', plot_fieldpairs{1}, plot_fieldpairs{2});
    exportgraphics(gcf, fullfile(fig_path,fig_name));
    
end