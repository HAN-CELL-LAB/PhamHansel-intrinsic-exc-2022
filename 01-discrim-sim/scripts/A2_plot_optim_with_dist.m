clc; clear; close all; 
run startup.m

%% Define path 
data_path = 'data/varyNYk'; 
fig_path = 'figures'; 

data_file = fullfile(data_path, 'summary_analyses.mat');
load(data_file, 'config', 'summary_analyses'); 

%% Get out variables and results
N_X = config.N_X;
N_Y_vec = config.N_Y_vec;
k_vec = config.k_vec;
dtheta_vec = config.dtheta_vec;

num_NY = length(N_Y_vec); 
num_k = length(k_vec); 

%% Plot
plt_select = struct(...
    'Jd', 'Jd_mean_nonan', ...
    'Hd', 'Hd_mean', ...
    'Y', 'Y_mean');

plt_fields = fieldnames(plt_select);

cmap = struct(...
    'Jd', return_colorbrewer('Reds', num_k) * 0.9, ...
    'Hd', return_colorbrewer('Blues', num_k) * 0.9, ...
    'Y', return_colorbrewer('Greys', num_k) * 0.9);

lbl_names = struct(...    
    'Jd', 'Jaccard', ...
    'Hd', 'Hamming', ...
    'Y', 'output');

graphic_setdefault(15, ...
    'DefaultAxesMinorGridAlpha', 0.05, ...
    'DefaultAxesMinorGridLineStyle', '-', ...
    'DefaultTextInterpreter', 'latex', ...
    'DefaultLegendInterpreter', 'latex', ...
    'DefaultStemMarkerSize', 1, ...
    'DefaultStemlineWidth', 1.5, ...
    'DefaultFigureWindowStyle','normal');


figure;

plt_netfields = {'init', 'best'};
for i_nf = 1:length(plt_netfields)
    
plt_netf = plt_netfields{i_nf}; 
anly_sel = summary_analyses.(plt_netf);
for i_ny = 1:num_NY
    subplot(2,3,i_ny + 3*(i_nf-1)); hold on;

    for i_k = 1:num_k
        k = k_vec(i_k);
        yyaxis left;
        plot(dtheta_vec, squeeze(anly_sel.(plt_select.Jd)(i_ny,i_k,:)), '-', 'linewidth', 3, ...
            'color', [cmap.Jd(i_k,:), 0.9])
        plot(dtheta_vec, squeeze(anly_sel.(plt_select.Hd)(i_ny,i_k,:)), '-',  'linewidth', 3, ...
            'color', [cmap.Hd(i_k,:), 0.9])
        
        yyaxis right;
        plot(dtheta_vec, squeeze(anly_sel.(plt_select.Y)(i_ny,i_k,:)), '-', 'linewidth', 2, ...
            'color', [cmap.Y(i_k,:), 0.9])        
        
    end
    
    yyaxis left; 
    set(gca, 'ycolor', 'k')
    if i_ny == 1
        ylabel('$J_d, H_d$ (colors)')
    end
        
    
    yyaxis right;
    set(gca, 'ycolor', 'k')
    
    if i_ny == num_NY
        ylabel('$Y$ (black)');
    end
    ylim([0,1])
    
    xlabel('$\Delta \theta$');
    title(sprintf('$N_Y = %d$ (%s)',  N_Y_vec(i_ny), plt_netf));
    
    yyaxis left;
end

end

linkaxes(findall(gcf, 'type', 'axes'), 'xy')

for i_f = 1:length(plt_fields)
    subplot(2,3,i_f); hold on;

    fn = plt_fields{i_f};
    colormap(gca, cmap.(fn));
    cbar = colorbar(gca, 'Location', 'west', ...
        'Ticks', linspace(0,1,num_k), 'TickLabels', k_vec);
    xlabel(cbar, lbl_names.(fn));
    title(cbar,'$k$', 'Interpreter', 'latex');
    cbar.Position = cbar.Position .* [1,1,0.5,0.5] + [0.005,0.05,0,0];

end

exportgraphics(gcf, fullfile(fig_path, 'dist-plot.pdf'));
