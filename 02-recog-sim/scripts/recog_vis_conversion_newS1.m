run start_up.m;
%% Load data
load('data/discrim_vs_recog_data.mat', ...
    'results_all', 'results_dim_description', ...
    'def_opts', 'param_opts');

fig_path = fullfile('figures/discrim-vs-recog/official');
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

%% Get range of Vthres
Vthres_range = [...
    mean([gill_s1_data.pooled_table.Vrest_1_base;gill_s1_data.pooled_table.Vrest_1_post]), ...
    mean(extra_s1_data.Vthres_max)];
 
%% Plot lines
common_norm_dthres_range = [-0.5, 0.5]; 

plot_fieldpairs = {'recognition_postradeoffs', 'recognition_TPR'};
title_name = struct; 
title_name.recognition_postradeoffs = 'recog. TPR - FPR';
title_name.recognition_TPR = 'recog. TPR';

cbar_name = struct; 
cbar_name.recognition_postradeoffs = 'recog. TPR - FPR';
cbar_name.recognition_TPR = 'recog. TPR';

cmap_pairs = {return_colorbrewer('Reds', length(percent_complete_input_vec)) * 0.9, ...
    return_colorbrewer('Blues', length(percent_complete_input_vec)) * 0.9};

nearest_select_alpha = [0.6, 0.75, 0.9];
nearest_select_novrl = [0];

select_alpha_ind = arrayfun(@(x) find_nearest(alpha_W_vec, x, 'ind'), nearest_select_alpha);
select_novrl_ind = arrayfun(@(x) find_nearest(num_overlap_per_group_vec, x, 'ind'), nearest_select_novrl);

select_alpha_vec = alpha_W_vec(select_alpha_ind);
select_novrl_vec = num_overlap_per_group_vec(select_novrl_ind);

fig_ncols = length(select_alpha_ind);
fig_nrows = length(select_novrl_vec);

plot_resultlist = {results_all.(plot_fieldpairs{1}), ...
    results_all.(plot_fieldpairs{2})};


figure('units', 'normalized', ...
    'Position', [0.05,0.05,1,0.7], ...
    'DefaultAxesFontSize', 20);

cnt_splt = 1;
pseudo_xaxis_pos = cell(length(select_novrl_ind),1);

fig_description = sprintf('\\color[rgb]{%f, %f, %f}%s\\color{black} vs \\color[rgb]{%f, %f, %f}%s', ...
    cmap_pairs{1}(end,:), title_name.(plot_fieldpairs{1}), ...
    cmap_pairs{2}(end,:), title_name.(plot_fieldpairs{2}));

annotation('textbox', 'String', fig_description, ...
    'FontSize', 30, 'Interpreter', 'tex', 'LineStyle', 'none', ...
    'Position', [0, 0.92, 1, 0.2], ...
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
        title(sprintf('$\\alpha_W = %.2f$', select_alpha_vec(i)));
        if j == 1
            pseudo_xaxis_pos{i} = get(gca, 'Position'); 
        end
        
        xlim(common_norm_dthres_range);
        despline;
        
        cnt_splt = cnt_splt + 1;
    end
    
    
end


for ax = findobj(gcf, 'type', 'axes')'
    ax.Position = ax.Position .* [1,1,1,0.5] + [0,0.3,0,0];
end

last_ax_pos = get(ax, 'Position');

common_text_opts_pseudo_xaxis = { ...
    'Fontsize', 20, 'Interpreter', 'latex',...
    'HorizontalAlignment', 'center'};

txtfmt_base = {'base %.1f mV', '%.1f mV'};
txtfmt_rest = {'%.1f mV', '%.1f mV'};
txtfmt_meanchange = {'%.2f mV', '%.2f mV'}; 
txtfmt_maxchange = {'%.2f mV', '%.2f mV'}; 

chosen_txtfmt = @(txtfmt_list, i) txtfmt_list{ (i~=1) + 1 };
changed_color = [0.1,0.4,0.2]; 
ticklength_scale = 0.25; 
txt_rng_YCoord = 0.5;
ax_YLim = [-1,5];

for i = 1:length(pseudo_xaxis_pos)
    ax_pos = pseudo_xaxis_pos{i};
    ax_pos = ax_pos .* [1,1,1,0.6] + [0,0.02,0,0];
    pseudo_xaxis = axes(gcf, 'Units', 'normalized', ...
        'Position', ax_pos, 'Visible', 'off');
    hold(pseudo_xaxis, 'on');

    pseudo_xaxis = gca;
    hold(pseudo_xaxis,'on');
    
    tickdat_x = common_norm_dthres_range;
    tickdat_y = zeros(1,length(tickdat_x));
    txtfmt_base_ith = chosen_txtfmt(txtfmt_base,i);
    txtfmt_rest_ith = chosen_txtfmt(txtfmt_rest,i);

    plot(pseudo_xaxis, tickdat_x([1,end]), tickdat_y([1,end]), ...
        'Color', 'k', 'LineWidth', 1.5);
    plot(pseudo_xaxis, [tickdat_x;tickdat_x], tickdat_y + ticklength_scale*[-0.5;0.5], ...
        'Color', 'k', 'LineWidth', 2);
    
    arrayfun(@(x,s,fmt) text(pseudo_xaxis, x, txt_rng_YCoord, sprintf(fmt{:}, s),...
        'Color', 'k', common_text_opts_pseudo_xaxis{:}), ...
        [-0.5,0.5], Vthres_range, ...
        {txtfmt_rest_ith, txtfmt_rest_ith});
      
    
    xlim(pseudo_xaxis, common_norm_dthres_range+[-0.02,0.02]);
    ylim(ax_YLim);
    despline;
end

linkaxes(findobj(gcf, 'type', 'axes'), 'x');
exportgraphics(gcf, fullfile(fig_path,'test_new_s1.pdf'));


