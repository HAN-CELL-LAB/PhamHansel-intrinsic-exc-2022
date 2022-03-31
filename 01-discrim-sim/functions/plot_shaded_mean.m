function [hobj_main, hobj_shade] = plot_shaded_mean(x, y_struct, z_sem, ...
                                        main_color, fill_alpha, additional_plotstyles, ...
                                        legend_label, sm_win, ...
                                        show_max, maxest_win)
    if ~exist('main_color', 'var'), main_color = [0.2,0.2,0.4]; end
    if ~exist('fill_alpha', 'var'), fill_alpha = 0.1; end
    if ~exist('additional_plotstyles', 'var'), additional_plotstyles = '-'; end
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
    hobj_main = plot(x, y_mean, 'color', main_color, additional_plotstyles{:});
    
    if ~isempty(legend_label) 
        set(hobj_main, 'DisplayName', legend_label, 'Tag', 'showlegend');
    end

end
