clc; clear; close all; 
run start_up; 

fig_path = fullfile('figures/toy-recog');
if ~exist(fig_path, 'dir')
    mkdir(fig_path);
end

%% Precursors for Fig1's
N = 10;
num_thres = 100;

%% Plot examples
graphic_reset(25, ...
    'DefaultAxesMinorGridAlpha', 0.05, ...
    'DefaultAxesMinorGridLineStyle', '-', ...
    'DefaultTextInterpreter', 'latex', ...
    'DefaultLegendInterpreter', 'latex', ...
    'DefaultStemMarkerSize', 1, ...
    'DefaultAxesTitleFontSize', 1.3, ...
    'DefaultAxesLabelFontSize', 1.2, ...
    'DefaultStemlineWidth', 1.5);

willplot = true;

plot_examples_combo_at_X(N, num_thres, 'example-equal', struct('type', 'equal'), willplot, fig_path, false);
plot_examples_combo_at_X(N, num_thres, 'example-uniform', struct('type', 'uniform'), willplot,fig_path, false);

%% Plot summary
close all; 
num_examples = 50;
equal_results = plot_examples_combo_at_X(N, num_thres, 1, struct('type', 'equal'));
uniform_results = arrayfun(@(id) plot_examples_combo_at_X(N, num_thres, id, struct('type', 'uniform')), 1:num_examples);


plot_names = fieldnames(equal_results);
plot_names = plot_names(1:2);
x_vecs = struct('W', 1:N, 'nYact', linspace(0,1,num_thres), 'Ypre', 0:2^N-1);
axes_labels = struct();
axes_labels.xlabel = struct('W', 'sorted input weight', 'nYact', 'threshold $\theta$',              'Ypre', 'combination ids');
axes_labels.ylabel = struct('W', '$W$',                 'nYact', '\# combinations',                 'Ypre', '$\sum_i W_i x_i$');
axes_labels.title  = struct('W', 'sorted weights',      'nYact', '\# of combinations with $Y=1$',   'Ypre', 'sorted weight sum');
unf_color = [0.2,0.2,0.2];

nrows = length(plot_names);

graphic_reset(22, ...
    'DefaultAxesMinorGridAlpha', 0.05, ...
    'DefaultAxesMinorGridLineStyle', '-', ...
    'DefaultTextInterpreter', 'latex', ...
    'DefaultLegendInterpreter', 'latex', ...
    'DefaultStemMarkerSize', 1, ...
    'DefaultAxesTitleFontSize', 1.25, ...
    'DefaultAxesLabelFontSize', 1.1, ...
    'DefaultStemlineWidth', 1.5);


shift_vertaxpos = [0.02; -0.02];
figure('units', 'normalized', 'position', [0,0,0.4,1]);
for i = 1:nrows
    subplot(nrows,1,i); hold on;
    
    field_name = plot_names{i};
    x_vec = x_vecs.(field_name);
    
    lgnd_objs = gobjects(2,1);
    lgnd_objs(1) = plot(x_vec, equal_results.(field_name), ':k', 'displayname', 'equals');
    
    unfres = vertcat(uniform_results.(field_name));
    mean_unfres = mean(unfres, 1);
    sem_unfres = std(unfres, 0, 1) / sqrt(num_examples);
    lower_bound = mean_unfres - sem_unfres;
    upper_bound = mean_unfres + sem_unfres;
    fill([x_vec fliplr(x_vec)], [lower_bound fliplr(upper_bound)], unf_color,...
        'LineStyle', 'none', 'FaceAlpha', 0.1)
    lgnd_objs(2) = plot(x_vec, mean_unfres, '--', 'color', unf_color, 'displayname', 'uniform');
    
    
    legend(lgnd_objs);
    
    ax_pos = get(gca, 'position');
    ax_pos(1) = ax_pos(1) + 0.05;
    ax_pos(2) = ax_pos(2) + shift_vertaxpos(i);
    set(gca, 'position', ax_pos);
    cellfun(@(lbl_id) set(get(gca, lbl_id), 'string', axes_labels.(lbl_id).(field_name)), fieldnames(axes_labels));
    
end

despline('all')

fig_name = fullfile(fig_path, 'summary.pdf');
exportgraphics(gcf, fig_name);
pause(0.5); close;
