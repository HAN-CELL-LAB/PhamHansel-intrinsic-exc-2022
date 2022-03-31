
function plot_single_result_panel(N_Y_vec, delta_theta_Y_vec, results, result_field, ...
    z_sem, cmap, fill_alpha, additional_plotstyles, ttl, ylbl, showlegend, lgnd_suff, lgnd_opts)
    if ~exist('showlegend', 'var'), showlegend = true; end
    if ~exist('lgnd_suff', 'var'), lgnd_suff = ''; end
    if ~exist('ylbl', 'var'), ylbl = result_field; end
    if ~exist('lgnd_opts', 'var'), lgnd_opts = {}; end 
    
    num_N_Y = length(N_Y_vec); 
    
    hold on; 
    for i = num_N_Y:-1:1
        main_color = cmap(i,:); 
        lgdn_lbl = sprintf('%d%s', N_Y_vec(i), lgnd_suff); 
        plot_shaded_mean(delta_theta_Y_vec, results{i}.(result_field), ...
            z_sem, main_color, fill_alpha, additional_plotstyles, lgdn_lbl);
    end

    xlabel('$\Delta \theta_Y$');
    ylabel(ylbl);
    title(ttl);
    if showlegend
        legend(findobj(gca, 'Tag', 'showlegend'), lgnd_opts{:})
    end
end