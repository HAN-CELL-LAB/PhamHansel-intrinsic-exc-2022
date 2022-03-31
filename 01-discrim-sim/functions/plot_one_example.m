function plot_one_example(results, prog_colors, prog_linestyles)
    if ~exist('prog_colors', 'var')
        prog_colors = struct('init', [0.75,0.75,0.95], 'best', [0.15,0.15,0.9]);
    end
    
    if ~exist('prog_linestyles', 'var')
        prog_linestyles = struct('init', '-', 'best', '-');
    end
    prog_fields = fieldnames(prog_colors);

    N_X = results.sim_config.N_X;
    N_Y = results.sim_config.N_Y;
    dtheta_vec = results.sim_config.dtheta_vec;
    
    figure;
    % optim prog
    subplot(221); hold on;
    plot(results.optim_progress.loss, '-k', 'displayname', 'loss', 'tag', 'showlegend'); 
    ylabel('loss');

    yyaxis right; set(gca, 'ycolor', 'k'); 
    plot(results.optim_progress.entropy, '-', 'color', [0.7,0.7,0.7], 'displayname', 'entropy', 'tag', 'showlegend');
    ylabel('entropy');
    xlim([-100, Inf]);

    cellfun(@(k) xline(results.net_config.(k).iter, '--', 'linewidth', 2, 'color', prog_colors.(k)), prog_fields);
    cellfun(@(k) plot(results.net_config.(k).iter, results.net_config.(k).entropy, ...
        'o', 'markersize', 8, 'color', prog_colors.(k), 'markerfacecolor', prog_colors.(k)), prog_fields);
    xlabel('iter');
    title(sprintf('optim progress ($%d X, %d Y$)', results.sim_config.N_X, results.sim_config.N_Y));
    legend(findobj(gca, 'tag', 'showlegend'));

    % compare init best net
    subplot(6,2,[7,9]); hold on;
    W0 = results.net_config.init.W;
    Wb = results.net_config.best.W;
    max_absW = max(abs([W0(:);Wb(:)]), [], 'all');

    image(1:N_X, 1:N_Y, W0);
    image((1:N_X) + N_X+2, 1:N_Y, Wb);
    colormap(return_colorbrewer('RdBu_R', 50));
    ax = gca; 
    ax.Position = ax.Position + [-0.03,0,0,0.05];
    cbar = colorbar('eastoutside');
    cbar.Position = cbar.Position .* [1,1,0.5,0.5] + [0.1,0.05,0,0];
    caxis([-1,1]*max_absW);
    daspect([1,1,1]);

    set(gca, 'TickLabelInterpreter', 'latex', 'XTick', [N_X/2, 2+N_X*3/2], ...
        'XTickLabel', {'$W_{YX}^{\mathrm{init}}$', '$W_{YX}^{\mathrm{best}}$'})
    despline();
    hide_only_axis('y');

    subplot(6,2,11); hold on;
    cellfun(@(k) plot(results.net_config.(k).theta, '.-', 'linewidth', 2, 'markersize', 30, ...
        'color', prog_colors.(k)), prog_fields);
    xlabel('$Y$ units'); ylabel('$\theta$')

    % compare ent, unq
    subplot(222); hold on; 
    cellfun(@(k) plot(dtheta_vec, results.apply_dtheta.(k).entropy, 'linewidth', 3, ...
        'color', prog_colors.(k), 'linestyle', prog_linestyles.(k), 'displayname', k), prog_fields);
    title('entropy '); xlabel('$\Delta \theta$'); ylabel('bits')
    legend('show'); despline();

    subplot(224); hold on; 
    cellfun(@(k) plot(dtheta_vec, results.apply_dtheta.(k).num_unq, 'linewidth', 3, ...
        'color', prog_colors.(k), 'linestyle', prog_linestyles.(k), 'displayname', k), prog_fields);
    title('num unq'); xlabel('$\Delta \theta$'); ylabel('\# unq patterns')
    legend('show'); despline();

end
