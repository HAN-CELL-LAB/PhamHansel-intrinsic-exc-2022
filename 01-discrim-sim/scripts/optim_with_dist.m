clc; clear; close all; 

addpath(genpath('extpkg'));
addpath(genpath('functions'));
graphic_setdefault(15, ...
    'DefaultAxesMinorGridAlpha', 0.05, ...
    'DefaultAxesMinorGridLineStyle', '-', ...
    'DefaultTextInterpreter', 'latex', ...
    'DefaultLegendInterpreter', 'latex', ...
    'DefaultStemMarkerSize', 1, ...
    'DefaultStemlineWidth', 1.5, ...
    'DefaultFigureWindowStyle','normal');

%%
N_X = 10;
N_Y_vec = [6, 10, 14]; 
k_vec = 2:8;

n_dtheta = 100;
n_sim = 100; 
n_iter = 2500; 

num_NY = length(N_Y_vec);
num_k = length(k_vec);
dtheta_vec = linspace(-1, 1, n_dtheta); 

Xids = 0:(2^N_X-1);
Xs = de2bi(Xids)';
  
% analyses = cell(num_NY, num_k, n_dtheta, n_sim);
% analyses(:) = {struct()}; 
% 
% % each N_Y
% for i_ny = 1:num_NY
%     N_Y = N_Y_vec(i_ny);
%     tic;
%     
%     ppm = ParforProgMon(sprintf('N_Y=%d', N_Y), n_sim);
%     % run multiple sims
%     parfor i_sim = 1:n_sim
%         
%         % optim
%         net_config = run_1_optim(N_X, N_Y, n_iter);
%         net_fields = {'init', 'best'};
%         
%         % either init or best 
%         for i_net = 1:length(net_fields)
%             nf = net_fields{i_net};
%             
%             % get net config
%             W = net_config.(nf).W;
%             theta = net_config.(nf).theta;
%             
%             % apply dtheta
%             for i_dtheta = 1:n_dtheta
%                 d_theta = dtheta_vec(i_dtheta);
%                 theta_app = theta + d_theta;
%                 
%                 Ys = W * Xs > theta_app;
%                 
%                 for i_k = 1:num_k
%                     k = k_vec(i_k);
%                     sub_X = sum(Xs(1:k,:), 1) == k;
%                     Y_subX = Ys(:,sub_X)';
%                     
%                     Jd = get_binary_pairwise_measure_distribution(Y_subX, @calculate_Jaccard_distance);
%                     Hd = get_binary_pairwise_measure_distribution(Y_subX, @calculate_Hamming_distance)/N_Y;
%                     
%                     tmp = struct(...
%                         'k', k, ...
%                         'Y_mean', mean(Y_subX,'all'), ...
%                         'Jd_mean', mean(Jd), ...
%                         'Jd_med', median(Jd), ...
%                         'Jd_sd', std(Jd), ...
%                         'Hd_mean', mean(Hd), ...
%                         'Hd_med', median(Hd), ...
%                         'Hd_sd', std(Hd),...
%                         'Jd_mean_nonan', mean(Jd, 'omitnan'), ...
%                         'Jd_med_nonan', median(Jd, 'omitnan'), ...
%                         'Jd_sd_nonan', std(Jd, 'omitnan') ...
%                         );
%                     
%                     analyses{i_ny, i_k, i_dtheta, i_sim}.(nf) = tmp; 
%                 end
%                 
%                 
%             end
%         end
%         
%         ppm.increment();
%     end
%     toc;
%     delete(ppm);
% end

%% 
% results = structarray_to_struct(cell2mat(results));
% config = struct(...
%     'N_X', N_X, ...
%     'N_Y_vec', N_Y_vec, ...
%     'k_vec', k_vec, ...
%     'n_dtheta', n_dtheta, ...
%     'n_sim', n_sim, ...
%     'n_iter', n_iter, ...
%     'dtheta_vec', dtheta_vec);

% save('data/varyNYk/sim_results.mat', 'analyses', 'config');
%%
% results = cell2mat(analyses);
% results = structarray_to_struct(results);
% results = structfun(@structarray_to_struct, results, 'uni', 0);

% anly = structfun(@(s) structfun(@(x) mean(x,4), s, 'uni', 0), results, 'uni', 0);

%%
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
% figure('units', 'normalized', 'position', [0, 0.2, 1, 0.5]);


plt_netfields = {'init', 'best'};
for i_nf = 1:length(plt_netfields)
    
plt_netf = plt_netfields{i_nf}; 
anly_sel = anly.(plt_netf);
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

exportgraphics(gcf, fullfile('figures/varyNYk', 'dist-Jdnonan.pdf'), 'ContentType', 'vector')