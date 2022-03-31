classdef multiClassifyObjA < classifyObjA
    properties
        num_weights
        original_props
    end
    methods
        function obj = multiClassifyObjA(props, num_weights)
            obj = obj@classifyObjA(props);
            obj.original_props = props; 
            obj.num_weights = num_weights;
        end
        
        function run(obj)
            stats_aggregrate = cell(obj.num_weights,1);
            props = obj.original_props; 
            for i = 1:obj.num_weights
                tmp_obj = classifyObjA(props);
                tmp_obj.run;
                stats_aggregrate{i} = tmp_obj.stats;
            end
            stats_aggregrate = vertcat(stats_aggregrate{:});
            
            obj.stats.fields = stats_aggregrate(1).fields; 
            for name_stats = ["all", "basedon_numactiveX"]
                stat_of_interest = [stats_aggregrate.(name_stats)]; 
                cat_stat = cellfun(@(name) cat(3,stat_of_interest.(name)), obj.stats.fields, 'uni', 0); 
                cat_stat = cell2struct(cat_stat, obj.stats.fields);
                obj.stats.(name_stats) = structfun(@(x) struct(...
                    'mean', mean(x,3), ...
                    'std', std(x,0,3)), cat_stat, 'uni', 0);
            end
            % save the last one as example
            obj.W_YX = tmp_obj.W_YX;
            obj.example_dataset = tmp_obj.example_dataset; 
        end
        
        function plot(obj, save_filename, save_in_par)
            if ~exist('save_in_par', 'var') 
                save_in_par = false;
            end
            %% Values to plot
            N_X = obj.N_X;
            thres_vec = obj.thres_vec;
            stats_all = obj.stats.all;
            
            stats_basedon_numactveX = obj.stats.basedon_numactiveX;
            num_unqactiveX = obj.num_unqactiveX;
            num_thres = obj.num_thres;
            
            W_YX = obj.W_YX;
            sz_objA = [1,1]*obj.N_side;
            W_YX_square = reshape(W_YX, sz_objA);
            [W_YX_sorted, ind_W_sort] = sort(W_YX, 'descend');
            ind_W_of_objA = obj.X_A == 1;
            ind_W_of_objA = ind_W_of_objA(ind_W_sort);
            range_W = [0, ceil(1e2*W_YX_sorted(1))/1e2];
            examples = obj.example_dataset;
            
            %% Plot setup 
           
            
            if save_in_par
                figure('units', 'inches', 'position', [0,0,17.708,9.979]);
            else
                figure('units', 'normalized', 'position', [0,0,0.88,0.9]);
            end
            
            plot_stats = {'TPR', 'FPR', 'TNR', 'FNR'};
            cmap_input = 0.95*gray(100);
            cmap_stats = return_colorbrewer('GnBu', 100);
            colororder_stats = return_colorbrewer('Paired');
            cmap_roccurve = return_colorbrewer('Spectral', num_thres)*0.98;
            
            nrows_splt = 4;
            ncols_splt = 6;
            
            input_image_plots = {...
                struct('subplot', 2, 'val', examples.objA_full, 'title', 'obj-\textbf{A}'), ...
                struct('subplot', 2+ncols_splt, 'val', examples.objA_input, 'title', 'example input'), ...
                struct('subplot', 3, 'val', examples.notA_full, 'title', 'not-\textbf{A}'), ...
                struct('subplot', 3+ncols_splt, 'val', examples.notA_input, 'title', 'example input'), ...
                };
            
            stats_image_spltidx = [5,6,5,6]+[0,0,1,1]*ncols_splt;

            %% 0. Annotation
            pseudo_ax = axes('Position', [0.14,0.6,0.1,0.1]);
            pseudo_ax.Visible = 'off'; 
            xlim(pseudo_ax,[0,1]); ylim(pseudo_ax,[0,1]);             
            
            ttl_description = '\textbf{DESCRIPTION}';
            netnum_description = sprintf(['$N_X = %d \\ (%d \\times %d)$' newline '$N_Y = %d$'], ...
                obj.N_X, obj.N_side, obj.N_side, obj.N_Y);
            
            objgen_description = sprintf(['obj\\textbf{A}: \\texttt{%s}', newline, ...
                'not\\textbf{A}: \\texttt{%s}', newline, 'dataset-size = %d'], ...
                lower(obj.name_of_objA),  lower(obj.name_of_notA), obj.size_input_dataset); 
            
            weight_description = sprintf(['gen-weight-meth: \\texttt{%s}', newline, ...
                'from-obj\\textbf{A}: $f=%s$', newline, ...
                'from-not\\textbf{A}: $f=%s$, $\\alpha=%.2f$', newline, ...
                'num-rnd-weights = %d'], ...
                obj.genWmethod, obj.genWmain_name, obj.genWnoise_name, obj.alpha_W, obj.num_weights);

            fig_description = strjoin({ttl_description, netnum_description, objgen_description, weight_description}, newline); 
            text(pseudo_ax, 0, 1, fig_description, 'FontSize', 10, ...
                'HorizontalAlignment', 'left', 'VerticalAlignment', 'cap');
            
            
            %% 1. Plotting examples objects and inputs             
            for i = 1:length(input_image_plots)
                plot_struct = input_image_plots{i};
                subplot(nrows_splt, ncols_splt, plot_struct.subplot); hold on;
                image_with_strict_limits(plot_struct.val);
                
                set(gca, 'xtick', '', 'ytick', '', 'box', 'on');
                colormap(gca, cmap_input); daspect([1,1,1]);
                title(plot_struct.title);
                
                if i == 2
                    xlabel('gray $\sim$ on');
                end
            end
            
            %% 2. Plotting example weight matrix
            subplot(nrows_splt, ncols_splt, 4); hold on;
            image_with_strict_limits(W_YX_square);
            set(gca, 'xtick', '', 'ytick', '', 'box', 'on');
            colormap(gca, 'gray'); daspect([1,1,1]); caxis(range_W);
            title('$W_{Y \to X}$');
            xlabel(sprintf('$N_X = %d$', N_X));
            ylabel(sprintf('$(%d \\times %d)$', sz_objA(1),sz_objA(2)));
            
            cbar = colorbar; ax_pos = get(gca, 'Position');
            cbar.Position = cbar.Position.*[1,0,0.8,0.7] + [0.03,ax_pos(2),0,0];
            cbar.Ticks = range_W;
            
            subplot(nrows_splt, ncols_splt, 4+ncols_splt); hold on;
            stem(find(~ind_W_of_objA), W_YX_sorted(~ind_W_of_objA), 'Color', [0.1,0.1,0.1]);
            stem(find(ind_W_of_objA), W_YX_sorted(ind_W_of_objA), ...
                'color', [0.8,0.8,0.8]);
            title('sorted weight'); xlabel('(gray is from A)')
            hide_only_axis('xy');
            
            %% 3. Plot mean "stats_basedon_numactveX"
            for i = 1:length(plot_stats)
                
                stat_to_plot = stats_basedon_numactveX.(plot_stats{i}).mean; 
                
                subplot(nrows_splt, ncols_splt, stats_image_spltidx(i)); hold on;
                image(thres_vec, num_unqactiveX, stat_to_plot);
                title(plot_stats{i});
                xlim(thres_vec([1,end])); ylim(num_unqactiveX([1,end]));
                caxis([0, 1]); colormap(gca, cmap_stats);
                pbaspect([1,1,1]);
                
                ax_pos = get(gca, 'Position'); 
                ax_pos(1) = ax_pos(1) + 0.03;
                set(gca, 'Position', ax_pos);
                
                if i ~= length(plot_stats)
                    set(gca, 'xtick', '', 'ytick', '', 'box', 'on');
                else
                    cbar = colorbar; ax_pos = get(gca, 'Position');
                    cbar.Position = cbar.Position.*[1,0,0.8,0.7] + [0.03,ax_pos(2),0,0];
                end
                if i == 3
                    xlabel('$\theta$ (threshold)');
                    ylabel('$k$ (\# active inputs)');
                end
            end
            
            %% 4. Plot "stats_all" varied by "threshold"
            axes('units', 'normalized', 'position', [0.2,0.1,0.4,0.3]); hold on;
            colororder(colororder_stats);
            lgnd_objs = cellfun(@(name) plot(thres_vec, stats_all.(name).mean, 'linewidth', 3, 'displayname', name), plot_stats);
            cellfun(@(name,color) fill([thres_vec, fliplr(thres_vec)], ...
                [stats_all.(name).mean - stats_all.(name).std, ...
                fliplr(stats_all.(name).mean + stats_all.(name).std)], ...
                color,  'LineStyle', 'none', 'FaceAlpha', 0.3), ...
                plot_stats, mat2colcell(colororder_stats(1:length(plot_stats),:))');
            ylim([0,1]);
            xlabel('$\theta$ (threshold)'); 
            ylabel(['mean (line) $\pm$ std (shade)' newline ...
                sprintf('(from %d random $$W_{YX}$)', obj.num_weights)]);
            title('performance rates', 'FontSize', 25);
            legend(lgnd_objs, 'Location', 'east');
            despline;
            
            %% 5. Plot "stats_all" ROC 
            axes('units', 'normalized', 'position', [0.7,0.1,0.2,0.3]); hold on;
            
            plot([0,1],[0,1], '--k');
            errorbar(stats_all.FPR.mean, stats_all.TPR.mean, ...
                 stats_all.TPR.std, stats_all.TPR.std, ...
                 stats_all.FPR.std, stats_all.FPR.std, ...
                 'Marker', 'none', 'Color', [0.9,0.9,0.9,0.1], 'LineWidth', 2, ...
                 'CapSize', 0.1);
            scatter(stats_all.FPR.mean, stats_all.TPR.mean, 30, cmap_roccurve,'filled', ...
                'MarkerFaceAlpha', 0.8, 'MarkerEdgeColor', 'none');
            xlabel('FPR'); ylabel('TPR'); title('ROC based on $\theta$');
            daspect([1,1,1]); despline;
            
            colormap(gca, cmap_roccurve);
            caxis(gca, thres_vec([1,end]));
            cbar = colorbar(gca);
            title(cbar, '\theta'); xlabel(cbar, 'threshold', 'Interpreter', 'latex');
            
            %% Save figure
            if exist('save_filename', 'var') 
                export_fig(gcf, save_filename, '-r300', '-p0.02', '-dpng');
            end
        end
    end
end