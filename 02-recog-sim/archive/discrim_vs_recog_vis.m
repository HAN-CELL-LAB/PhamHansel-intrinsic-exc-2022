run start_up.m;
%% Load data
load('data/discrim_vs_recog_data.mat', ...
    'results_all', 'results_dim_description', ...
    'def_opts', 'param_opts'); 

fig_path = fullfile('figures/discrim-vs-recog/misc');
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

%% Plot subplot heatmaps 
fig_nrows = length(percent_complete_input_vec);
fig_ncols = length(num_overlap_per_group_vec); 

cmap = return_colorbrewer('RdBu_r', 50) * 0.95; 

for ind_field = 1:length(plot_fields)
    plot_field = plot_fields{ind_field};
    plot_result = results_all.(plot_field);
    caxis_lim = [min(plot_result, [], 'all'), max(plot_result, [], 'all')]; 
    
    figure(...
        'units','normalized', ...
        'Position',[0.05,0.05,0.5,0.95], ...
        'DefaultAxesFontSize', 12);
    
    colormap(cmap);
    
    cur_splt = 1;
    annotation('textbox', 'String', ...
        sprintf('\\textbf{%s}', regexprep(plot_field, '_', '-')), ...
        'FontSize', 15, 'Interpreter', 'latex', 'LineStyle', 'none', ...
        'Position', [0.1, 0.92, 0.2, 0.05], ...
        'HorizontalAlignment', 'left', 'VerticalAlignment', 'cap');
    
    for ind_overlap = 1:length(num_overlap_per_group_vec)
        for ind_percent = 1:length(percent_complete_input_vec)
            
            percent_complete_input = percent_complete_input_vec(ind_percent);
            num_overlap_per_group = num_overlap_per_group_vec(ind_overlap);
            
            plot_mat = squeeze(plot_result(:, :, ind_percent, ind_overlap));
            
            subplot(fig_nrows, fig_ncols, cur_splt); hold on;
            image(dthres_vec, alpha_W_vec, plot_mat');

            xlim(dthres_vec([1,end])); ylim(alpha_W_vec([1,end]));
            set(gca, 'box', 'off');
            caxis(caxis_lim);
            pbaspect([1,1,1]);
            
            if ind_percent == 1
                ylabel({sprintf('$n_{ovrl} = %d$', num_overlap_per_group), '$\alpha_W$'});
            end
            
            if ind_overlap == length(num_overlap_per_group_vec)
                xlabel({'$\Delta\theta$', sprintf('$\\%%c = %.1f$', percent_complete_input)});
            end
            
            if ind_percent ~= 1 || ind_overlap ~= length(num_overlap_per_group_vec) 
                set(gca, 'xtick', '', 'ytick', '');
                hide_only_axis(gca, 'xy');
            end
            
            if ind_percent == length(percent_complete_input_vec) && ind_overlap == length(num_overlap_per_group_vec) 
                cbar = colorbar; 
                cbar.FontSize = 8;
                cbar.Position = cbar.Position .* [1,1,1.5,2.5] + [0.05,-0.02,0,0];
            end
            cur_splt = cur_splt + 1;
            
        end
    end
    
    fig_name = fullfile(fig_path, sprintf('%s-heatmapsubplot', plot_field));
    export_fig(fig_name,  '-r300', '-p0.02');
    pause(0.5); close; 
end

%% Plot 4D heatmaps
fig_nrows = 2; 
fig_ncols = 2;
padding_size = 2; 
select_plotfields = {'discrimination_accuracy', 'recognition_accuracy', ...
    'recognition_L2_truerates', 'recognition_L2_falserates'}; 

figure('units','normalized', 'Position',[0.05,0,0.7,1]);

for ind_field = 1:length(select_plotfields)
    plot_field = select_plotfields{ind_field};
    plot_result = results_all.(plot_field);
    plot_result = squeeze(mat2cell(plot_result, ...
        def_opts.num_dthres, ...
        length(alpha_W_vec), ...
        ones(1,length(percent_complete_input_vec)), ...
        ones(1,length(num_overlap_per_group_vec))));
    plot_result = cellfun(@(x) padarray(x, padding_size*[1,1], nan, 'post'), plot_result, 'uni', 0);
    plot_result = cell2mat(plot_result);
    
    subplot(fig_nrows, fig_ncols, ind_field); hold on;
    image(percent_complete_input_vec, num_overlap_per_group_vec, plot_result', ...
        'AlphaData', ~isnan(plot_result'));
    xlabel('\hspace{7em} \% complete');
    ylabel('\# overlap per group');
    title(regexprep(plot_field, '_', '-'));

    xlim(percent_complete_input_vec([1,end]));
    ylim(num_overlap_per_group_vec([1,end]));
    set(gca, 'xtick', percent_complete_input_vec, 'ytick', num_overlap_per_group_vec);
    
    colormap(cmap);
    pbaspect([1,1,1]);
    
    cbar = colorbar('Location', 'eastoutside', 'FontSize', 12);
    cbar.Position = cbar.Position  .* [1,1,0.5,0.5] + [0.01,0,0,0];
    despline;
    
end


mini_ax = axes('Position', [0.89,0.82,0.1,0.1], 'Color', [0.95,0.95,0.95]);
text(dthres_vec(1) + 1/20 * diff(dthres_vec([1,end])), ...
    alpha_W_vec(1) +  1/2 * diff(alpha_W_vec([1,end])), ...
    ['mini-axes' newline 'inside'], 'fontsize', 15);

xlim(mini_ax, dthres_vec([1,end]));
ylim(mini_ax, alpha_W_vec([1,end]));
xlabel('$\Delta\theta$'); 
ylabel('$\alpha_W$');
pbaspect([1,1,1]); 
despline(1.5);

  
fig_name = fullfile(fig_path, 'selectfields-heatmapflattend4D');
export_fig(fig_name,  '-r600', '-p0.02');
pause(0.5); close;

%% Plot lines 

% plot_fieldpairs = {'discrimination_accuracy', 'recognition_accuracy'}; 
% plot_fieldpairs = {'recognition_TPR', 'recognition_FPR'}; 
% plot_fieldpairs = {'recognition_L1_truerates', 'recognition_L1_falserates'}; 
% plot_fieldpairs = {'recognition_L2_truerates', 'recognition_L2_falserates'}; 
% plot_fieldpairs = {'recognition_TPR', 'recognition_postradeoffs'}; 
plot_fieldpairs = {'recognition_FPR', 'recognition_postradeoffs'}; 

% plot_fieldpairs = {'discrimination_accuracy', 'recognition_L2_truerates'};
% plot_fieldpairs = {'discrimination_accuracy', 'recognition_L2_falserates'};

% plot_fieldpairs = {'recognition_L1_truerates', 'recognition_L1_falserates'};
% plot_fieldpairs = {'recognition_TPR', 'recognition_FNR'};
% plot_fieldpairs = {'recognition_TNR', 'recognition_FPR'};

cmap_pairs = {return_colorbrewer('Reds', length(percent_complete_input_vec)) * 0.9, ...
    return_colorbrewer('Blues', length(percent_complete_input_vec)) * 0.9}; 

select_alpha_ind = [20,30,35,40,45,50]; 
select_alpha_vec = alpha_W_vec(select_alpha_ind); 

select_novrl_ind = [1,3,5,7]; 
select_novrl_vec = num_overlap_per_group_vec(select_novrl_ind); 

fig_ncols = length(select_alpha_ind); 
fig_nrows = length(select_novrl_vec);

plot_resultlist = {results_all.(plot_fieldpairs{1}), ...
    results_all.(plot_fieldpairs{2})}; 

figure; 
cnt_splt = 1;

fig_description = sprintf('\\color[rgb]{%f, %f, %f}%s\\color{black} vs \\color[rgb]{%f, %f, %f}%s', ...
    cmap_pairs{1}(end,:), regexprep(plot_fieldpairs{1}, '_', '-'), ...
    cmap_pairs{2}(end,:), regexprep(plot_fieldpairs{2}, '_', '-'));

annotation('textbox', 'String', fig_description, ...
    'FontSize', 20, 'Interpreter', 'tex', 'LineStyle', 'none', ...
    'Position', [0, 0.92, 1, 0.08], ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');

for j = 1:length(select_novrl_ind)
    for i = 1:length(select_alpha_ind)
        
        subplot(fig_nrows,fig_ncols,cnt_splt); hold on;
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

linkaxes(findall(gcf, 'type', 'axes'), 'xy'); 
% xlim([-0.75,0.75]); 
xlim([-0.5,0.5]); 

cbar_tickval = percent_complete_input_vec; 
cbar_numtick = length(cbar_tickval); 
cbar_tickpos = linspace(0.5/cbar_numtick,1-0.5/cbar_numtick,cbar_numtick);
cbar_props = {'ticks', cbar_tickpos, 'ticklabels', cbar_tickval};

last_ax_pos = get(gca, 'Position');
pseudo_ax_1 = axes('Position', last_ax_pos, 'Visible', 'off');
pseudo_ax_2 = axes('Position', last_ax_pos, 'Visible', 'off');

colormap(pseudo_ax_1, cmap_pairs{1});
cbar_1 = colorbar(pseudo_ax_1); 
cbar_1.Position = cbar_1.Position .* [1,1,1,0.5] + [0,0,0,0];
cbar_1.Position(1) = sum(last_ax_pos([1,3])) + 0.02; 
ylabel(cbar_1, ['% complete' newline sprintf('(%s)', regexprep(plot_fieldpairs{1}, '_', '-'))]); 
set(cbar_1, cbar_props{:});

colormap(pseudo_ax_2, cmap_pairs{2});
cbar_2 = colorbar(pseudo_ax_2); 
cbar_2.Position = cbar_1.Position + [0,0.25,0,0];
ylabel(cbar_2, ['% complete' newline sprintf('(%s)', regexprep(plot_fieldpairs{2}, '_', '-'))]); 
set(cbar_2,  cbar_props{:});

% fig_name = fullfile(fig_path, sprintf('%s-AND-%s-lineplotOfalphaWAndnOverlap-lambda1', plot_fieldpairs{1}, plot_fieldpairs{2}));
% export_fig(fig_name,  '-r300', '-p0.02');
% pause(0.5); close;


