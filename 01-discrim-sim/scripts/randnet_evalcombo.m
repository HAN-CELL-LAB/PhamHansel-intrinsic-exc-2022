clc; clear; close all; 

addpath(genpath('extpkg'));
graphic_setdefault(15, ...
    'DefaultAxesMinorGridAlpha', 0.05, ...
    'DefaultAxesMinorGridLineStyle', '-', ...
    'DefaultTextInterpreter', 'latex', ...
    'DefaultLegendInterpreter', 'latex', ...
    'DefaultStemMarkerSize', 1, ...
    'DefaultStemlineWidth', 1.5, ...
    'DefaultFigureWindowStyle','normal');


%%
if isempty(gcp('nocreate'))
    parpool('local', 6);
end

%%
N_X = 10; 
N_Y_vec = 1:20; 
theta_Y_base = 0.5; 
delta_theta_Y_vec = linspace(-0.5, 0.5, 200);
n_sim = 100; 

num_N_Y = length(N_Y_vec); 
results = cell(num_N_Y, 1); 
tic
parfor i = 1:num_N_Y
    N_Y = N_Y_vec(i); 
    results{i} = process_theta(N_X, N_Y, theta_Y_base, delta_theta_Y_vec, n_sim)
end
toc

%%
z_sem = 3; 
fill_alpha = 0.2;
cmap = jet(num_N_Y)*0.9; 

figure; 
subplot(211); hold on; 
ttl = sprintf('Entropy of output patterns (max = %d)', N_X);
plot_single_result_panel(N_Y_vec, delta_theta_Y_vec, results, 'entropy', ...
    z_sem, cmap, fill_alpha, ttl  , 'entropy (bits)', 1, {'NumColumns', 2'});

subplot(212); hold on; 
ttl = sprintf('Number of unique output patterns (max = %d)', 2^N_X);
plot_single_result_panel(N_Y_vec, delta_theta_Y_vec, results, 'num_unq', ...
    z_sem, cmap, fill_alpha, ttl, '\# patterns', false)
% set(gca, 'yscale', 'log')

%%
function plot_single_result_panel(N_Y_vec, delta_theta_Y_vec, results, result_field, ...
    z_sem, cmap, fill_alpha, ttl, ylbl, showlegend, lgnd_opts)
    if ~exist('showlegend', 'var'), showlegend = true; end
    if ~exist('ylbl', 'var'), ylbl = result_field; end
    if ~exist('lgnd_opts', 'var'), lgnd_opts = {}; end 
    
    num_N_Y = length(N_Y_vec); 
    
    hold on; 
    for i = num_N_Y:-1:1
        main_color = cmap(i,:); 
        lgdn_lbl = sprintf('%d', N_Y_vec(i)); 
        plot_shaded_mean(delta_theta_Y_vec, results{i}.(result_field), ...
            z_sem, main_color, fill_alpha, lgdn_lbl);
    end

    xlabel('$\Delta \theta_Y$');
    ylabel(ylbl);
    title(ttl);
    if showlegend
        legend(findobj(gca, 'Tag', 'showlegend'), lgnd_opts{:})
    end
end

function [hobj_main, hobj_shade] = plot_shaded_mean(x, y_struct, z_sem, ...
                                        main_color, fill_alpha, legend_label, sm_win, ...
                                        show_max, maxest_win)
    if ~exist('main_color', 'var'), main_color = [0.2,0.2,0.4]; end
    if ~exist('fill_alpha', 'var'), fill_alpha = 0.1; end
    if ~exist('legend_label', 'var'), legend_label = ''; end
    if ~exist('sm_win', 'var'), sm_win = 1; end
    if ~exist('show_max', 'var'), show_max = true; end
    if ~exist('maxest_win', 'var'), maxest_win = 50; end

    y_mean = smooth(y_struct.mean,sm_win);
    [~, pk_loc] = max(smooth(y_mean,maxest_win));

    y_err = z_sem * smooth(y_struct.sem,sm_win);
    y_low = y_mean - y_err;
    y_upp = y_mean + y_err;

    hold on; 
    hobj_shade = fill([x, fliplr(x)], [y_low; flipud(y_upp)], main_color, ...
        'FaceAlpha', fill_alpha, 'LineStyle', 'none'); 
    if show_max
        plot(x(pk_loc), y_mean(pk_loc), 'color', [main_color, 0.8], 'marker', '.', 'markersize', 30);
    end
    hobj_main = plot(x, y_mean, 'color', main_color);

%     hobj_main = plot(x(1:end-1), abs(diff(y_mean)), 'color', main_color);
%     hobj_main = plot(x, cumsum(y_mean), 'color', main_color);

    
    if ~isempty(legend_label) 
        set(hobj_main, 'DisplayName', legend_label, 'Tag', 'showlegend');
    end

end

function res = process_theta(N_X, N_Y, ...
                            theta_Y_base, delta_theta_Y_vec, ...
                            n_sim, take_stats)
    persistent Xs;
    if isempty(Xs)
        Xs = de2bi(0:(2^N_X-1))';
    end

    if ~exist('n_sim', 'var') 
        n_sim = 1; 
    end

    if ~exist('take_stats', 'var')
        take_stats = true;
    end

    n_dtheta = length(delta_theta_Y_vec); 

    res = struct(...
        'entropy', zeros(n_dtheta, n_sim), ...
        'num_unq', zeros(n_dtheta, n_sim) ...
        ); 

    for i_dtheta = 1:n_dtheta
        delta_theta_Y = delta_theta_Y_vec(i_dtheta); 
        theta_Y = config_theta(theta_Y_base, delta_theta_Y, N_Y);

        for i_sim = 1:n_sim
            W = rand(N_Y,N_X);
            W = W./sum(W,2);
            Ys = W * Xs > theta_Y;
            Yids = bi2de(Ys');

            [p,k] = id2prob(Yids);
            ent = entropy(p);

            res.entropy(i_dtheta, i_sim) = ent;
            res.num_unq(i_dtheta, i_sim) = k;
        end
    end

    if take_stats
        res = structfun(@(x) process_stats(x), res, 'uni', 0);
    end

end

function res = process_stats(X)
    n = size(X,2); 
    res = struct('mean', mean(X,2), 'sem', std(X,[],2)/sqrt(n));
end

function theta = config_theta(theta_base, delta_theta, N)
    theta = theta_base * ones(N, 1) + delta_theta; 
end

function [p, k] = id2prob(v) 
    p = groupcounts(v); 
    p = p / sum(p); 
    k = length(p); 
end

function ent = entropy(p)
    p = p / sum(p); 
    ent = -p .* log2(p); 
    ent = sum(ent(~isinf(ent) & ~isnan(ent))); 
end
