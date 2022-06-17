clc; clear; close all; 
run start_up; 

fig_path = fullfile('figures/demo');
if ~exist(fig_path, 'dir')
    mkdir(fig_path);
end

%% Activation function 
ip_color = [0.7,0.7,0.7,0.5]; 

x = -10:0.01:10; 
g = 1000; 
figure('position', [0,0,1,1]); hold on; 
plot(x, 1./(1+exp(-g*x)), '-k', 'linewidth', 10);
plot(x, 1./(1+exp(-g*(x - 5))), '-', 'color', ip_color,  'linewidth', 10);
plot(x, 1./(1+exp(-g*(x + 5))), '-', 'color', ip_color, 'linewidth', 10);

text(-3.5, 0.85, '$\uparrow$ IE', 'Color', ip_color, 'fontsize', 80, 'Interpreter', 'latex');
text(-3.5, 0.65, '$\downarrow \theta$', 'Color', ip_color, 'fontsize', 80, 'Interpreter', 'latex');

text(1, 0.2, '$\downarrow$ IE', 'Color', ip_color, 'fontsize', 80, 'Interpreter', 'latex');
text(1, 0.4, '$\uparrow \theta$', 'Color', ip_color, 'fontsize', 80, 'Interpreter', 'latex');

annotation('textarrow',[0.5, 1/3],[0.5, 0.5], ...
    'Color', ip_color, 'LineWidth', 5, ...
    'HeadLength', 50, 'HeadWidth', 50);

annotation('textarrow',[0.535, 0.7],[0.5, 0.5], ...
    'Color', ip_color, 'LineWidth', 5, ...
    'HeadLength', 50, 'HeadWidth', 50);


set(gca, 'xcolor', 'none', 'ycolor', 'none'); 
title('activation', 'fontsize', 100); 

exportgraphics(gcf, fullfile(fig_path, 'ip_heaviside.pdf'));

pause(1); close; 

%% Default options 
N_Y = 3;
N_X = N_Y * 50;
 
num_dthres = 10;
thres_baseline = 0.5;
num_inputs = 1;
sigma_noise = 0.8;
rng_dthres = [-1,1];

num_overlap_per_group = 0;
percent_complete_input = 1;
alpha_W = 0.5;

def_opts = struct(...
    'N_X', N_X, ...
    'N_Y', N_Y, ...
    'demo', true, ...
    'alpha_W', alpha_W, ...
    'num_dthres', num_dthres, ...
    'thres_baseline', thres_baseline, ...
    'num_inputs', num_inputs, ...
    'sigma_noise', sigma_noise, ...
    'rng_dthres', rng_dthres, ...
    'percent_complete_input', percent_complete_input, ...
    'num_overlap_per_group', num_overlap_per_group);

%% Essential words 
word_list = {...
    '$\Delta\theta$', '$\Delta\theta_{\mathrm{optim}}$', ...
    sprintf('$\\theta_{\\mathrm{base}}$ = %g', thres_baseline), ...
    sprintf('$\\sigma_{\\mathrm{noise}}$ = %g', sigma_noise), ...
    '$W_{YX}$', ...
    '$\alpha_W$', '\% complete', '\# overlap', ...
    'Example input (without noise)', ...
    'Expected output', ...
    'Discrimination', ...
    'Recognition'}; 

phrase_list = {...
    sprintf('Input depends on:\n(1) \\%% complete $\\in$ [0,1]\n(2) %s', ...
        sprintf('$\\sigma_{\\mathrm{noise}}$ = %g', sigma_noise)), ...
    sprintf('$W_{YX}$ depends on:\n(1) $\\alpha_W \\in [0,1] $\n(2) \\# overlap $$\\in [0,45]$$'), ...
    sprintf('Activation function depends on:\n(1) $\\Delta\\theta \\in [-0.5,0.5]$\n(2) %s', ...
        sprintf('$\\theta_{\\mathrm{base}}$ = %g', thres_baseline)), ...
    sprintf('Task:\n(1) Discrimination\n(2) Recognition'), ...
    ['Change of $\theta$ on output \#1' newline, ...
    '$\theta_1 = 0.5 + \Delta\theta$' newline, ...
        '$\theta_2 = \theta_3 = 0.5$ (base)']};

all_label_list = {word_list, phrase_list}; 

figure('position', [0,0,1,1]);

arrayfun(@(ind) ...
annotation('textbox', ...
    'Position', [0.05 + 0.5*(ind-1),0.05,0.45,0.9], ...
    'String', sprintf('%s\n', all_label_list{ind}{:}), ...
    'FontSize', 30, 'Interpreter', 'latex', ...
    'VerticalAlignment', 'middle', 'HorizontalAlignment', 'left', ...
    'LineStyle', 'none', 'FitBoxToText', 'on'), 1:length(all_label_list)); 


fig_name = fullfile(fig_path, 'words.pdf');
exportgraphics(gcf, fig_name);
pause(1); close;

%% Demonstrate ALPHA_W and N_OVERLAP
alpha_W_vec = [0.2, 0.4, 0.6, 0.8]; 
num_overlap_per_group_vec = [0, 20];

figure('Units', 'normalized', 'Position', [0.05,0.05,0.28,0.9]);
colormap('gray'); 
cnt_splt = 1; 
for i = 1:length(num_overlap_per_group_vec)
    for j = 1:length(alpha_W_vec)
        subplot(length(num_overlap_per_group_vec), length(alpha_W_vec), cnt_splt); 
        hold on; cnt_splt = cnt_splt + 1; 
        
        opts_tmp =  def_opts; 
        opts_tmp.num_overlap_per_group = num_overlap_per_group_vec(i); 
        opts_tmp.alpha_W = alpha_W_vec(j); 
        res_tmp = get_specperf(opts_tmp); 
        image_with_strict_limits(res_tmp.W_YX); 
        
        pbaspect([0.15,1,1]);
        
        ax_pos = get(gca, 'Position');
        ax_pos = ax_pos + [-0.01,0.05,0,0] * (i ~= 1);
        set(gca, 'position', ax_pos, 'box', 'on', 'linewidth', 2, 'xtick', '', 'ytick', '');
        caxis([0,0.02]);
        if j == 1
            ylabel([sprintf('\\# overlap = %d', opts_tmp.num_overlap_per_group) newline], 'FontSize', 18);
        end
        
        if i == length(num_overlap_per_group_vec)
            xlabel([newline sprintf('$\\alpha_W = %.2g$', opts_tmp.alpha_W)], ...
                'Rotation', 30,  'FontSize', 18, ...
                'VerticalAlignment', 'top', 'HorizontalAlignment', 'right');
        end
        
    end
end

last_ax_pos = get(gca, 'position'); 
cbar = colorbar('fontsize', 12, 'ticks', [0, 0.01, 0.02]);
cbar.Position(1) = last_ax_pos(1) + last_ax_pos(3) + 0.01;
cbar.Position(2) = last_ax_pos(2);
cbar.Position = cbar.Position .* [1,1,2.5,0.5];
subplot(length(num_overlap_per_group_vec), length(alpha_W_vec), 2); 
xlabel(sprintf('\\textbf{to} (%d Y)', N_Y), 'fontsize', 20);
ylabel(sprintf('\\textbf{from} (%d X)', N_X), 'fontsize', 20);

annotation('textbox', 'String', '$\mathbf{W_{YX}}$', 'fontsize', 35, ...
    'Units', 'normalized', 'Position', [0,0.94,1,0.05], ...
    'VerticalAlignment', 'middle', 'HorizontalAlignment', 'center', ...
    'LineStyle', 'none', 'Interpreter', 'latex')

fig_name = fullfile(fig_path, 'WXY_alphaW-noverlap.pdf');
exportgraphics(gcf, fig_name);
pause(1); close;

%% Demonstrate %COMPLETE
percent_complete_input_vec = [0.3,0.8,1]; 

output_colors = [200,40,40;40,160,40;40,40,200]/255;
normal_color = [0.7,0.7,0.7];

heading_label = '\color[rgb]{0,0,0} \rightarrow  '; 
spacing_label = '  '; 


nskip_splt_within_group = 1;
nskip_splt_between_group = 0;
nrows = 1; 
ncols = (sum(percent_complete_input_vec ~= 1) * def_opts.num_inputs + 1)* (N_Y + nskip_splt_within_group) + ...
    nskip_splt_between_group * length(percent_complete_input_vec); 
ttl_splt_order = [2,2,2]; 

init_recog_label   = repmat({'-'}, [1,N_Y]);
       
figure('Units', 'normalized', 'Position', [0.05,0.05,0.6,0.8]);
colormap('gray');

normal_color_str = sprintf('\\color[rgb]{%f,%f,%f}',normal_color);

cnt_splt = 1; 
for i = 1:length(percent_complete_input_vec)
    opts_tmp =  def_opts;
    opts_tmp.percent_complete_input = percent_complete_input_vec(i);
    if opts_tmp.percent_complete_input == 1
        opts_tmp.num_inputs = 1;
    end
    
    res_tmp = get_specperf(opts_tmp);
    
    masked_inputs = res_tmp.masked_inputs; 
    output_labels = res_tmp.labels; 
    for j = 1:length(masked_inputs)
        subplot(ncols,nrows,cnt_splt); hold on;
        if j == ttl_splt_order(i)
            ylabel(sprintf('\\%% complete = %g', opts_tmp.percent_complete_input), 'fontsize', 18);
        end
        image_with_strict_limits(masked_inputs{j});
        set(gca, 'xtick', '', 'ytick', '', 'box', 'on', 'linewidth', 3);
        
        cnt_splt = cnt_splt+1;
        if mod(j, N_Y) == 0 
            cnt_splt = cnt_splt + nskip_splt_within_group;
        end
        
        yyaxis right;
        ax_pos = get(gca, 'Position'); 
        ax_pos(3) = ax_pos(3) * 0.8;
        set(gca, 'position', ax_pos, 'ycolor', 'k', 'ytick', '');
        
        recog_label   = init_recog_label;
        
        highlighted_label = sprintf('\\color[rgb]{%f,%f,%f}1%s',...
            output_colors(output_labels(j),:), normal_color_str);
        
        recog_label{output_labels(j)}   = highlighted_label;
        
        if output_labels(j) ~= 1
            recog_label = {''};
        end
        
        label_combined = sprintf('%s%s%s%s%s', ...
            heading_label, ...
            normal_color_str, ...
            spacing_label, ...
            sprintf('%s',recog_label{:})); 
        ylabel(label_combined, 'FontSize', 30,  'Interpreter', 'tex', ...
            'Rotation',  0,  'VerticalAlignment', 'middle', 'HorizontalAlignment', 'left');
    end
    
    cnt_splt = cnt_splt + nskip_splt_between_group;
    
end

fig_name = fullfile(fig_path, 'inpout_percentcomplete.pdf');
exportgraphics(gcf, fig_name);
pause(1); close;