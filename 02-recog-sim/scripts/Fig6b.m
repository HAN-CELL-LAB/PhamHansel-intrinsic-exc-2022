clc; clear; close all; 
run start_up.m;

%% Load data
load('data/discrim_vs_recog_data.mat', ...
    'results_all', 'results_dim_description', ...
    'def_opts', 'param_opts');

fig_path = fullfile('figures/official');
if ~exist(fig_path, 'dir')
    mkdir(fig_path);
end

load('data/conversion_data_Li2020_individualcurved-morepoints.mat', ...
    'convert_data_obj', 'theta_range', 'theta_base', 'num_points');

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


%% Get optimal and covert

convert_obj = convert_data_obj.results.quadratic;
diff_from_optim_acc = 0.001;

select_result_fields = {'discrimination_accuracy', 'recognition_postradeoffs', 'recognition_TPR'};
optim_results = struct();
dtheta_optim = struct();
dVm_optim = struct; 
Vm_optim = struct;
for i = 1:length(select_result_fields)
    field_name = select_result_fields{i};
    result_curve = results_all.(field_name);
    dim_prm = size(result_curve);
    dim_prm = dim_prm(2:end); 
    
    perf_max = nan(dim_prm);
    dtheta_mat = nan(dim_prm);
    dVm_mat = nan(dim_prm);
    Vm_mat = nan(dim_prm); 
    
    
    for i1 = 1:dim_prm(1)
        for i2 = 1:dim_prm(2)
            for i3 = 1:dim_prm(3)
                Y_Vec = squeeze(result_curve(:,i1,i2,i3));
                
                diff_from_max_acc = abs(Y_Vec -  max(Y_Vec));
                dthres_discrim_max = dthres_vec(diff_from_max_acc < diff_from_optim_acc);
                dthres_optim = find_nearest(dthres_discrim_max, 0); 
                dtheta_mat(i1,i2,i3) = dthres_optim;
                perf_max(i1, i2, i3) = max(Y_Vec); 
                
                if abs(dthres_optim) <= 0.52
                    optim_ind_lookup = find_nearest(convert_obj.delta_theta_vs_base, dthres_optim, 'ind');
                    dVm_mat(i1,i2,i3) = convert_obj.delta_Vt_vs_base(optim_ind_lookup);
                    Vm_mat(i1,i2,i3) = convert_obj.Vt_vec(optim_ind_lookup);
                    
                end
            end
        end
    end
    
        
    optim_results.(field_name) = perf_max; 
    dtheta_optim.(field_name) = dtheta_mat;
    dVm_optim.(field_name) = dVm_mat;
    Vm_optim.(field_name) = Vm_mat;
end    


%% Select demo and plot
select_pairs = {...
    struct('alpha', 0.53, 'inp', 0.7), ...
    struct('alpha', 0.81, 'inp', 0.8), ...
    struct('alpha', 0.98, 'inp', 0.8)};

select_novlp_ind = 1; 
ind_percent_inp = 7:length(percent_complete_input_vec);

% smooth_fun = @(x) x;
smooth_fun = @(x) smoothdata(x, 'gaussian', 5);

num_colors = length(ind_percent_inp);
scale_color = 0.9;
cmaps = struct(...
    'discrim', return_colorbrewer('Reds', num_colors) * scale_color, ...
    'recognition_postradeoffs', return_colorbrewer('Greens', num_colors) * scale_color, ...
    'recognition_TPR', return_colorbrewer('Blues', num_colors) * scale_color);

plt_optim_order = {'dV', 'max'};

lgnd_shorts = struct(...
    'discrim', 'discrim$_{\mathrm{acc}}$', ...
    'recognition_postradeoffs', 'recog$_{\mathrm{TPR-FPR}}$', ...
    'recognition_TPR', 'recog$_{\mathrm{TPR}}$');

recogplt_opts = struct(...
    'recognition_postradeoffs', ...
        struct(...
            'title', 'discrim$_{\mathrm{acc}}$ and recog$_{\mathrm{TPR-FPR}}$', ...
            'demo_showlegend', true, ...
            'optim_showlegend', false, ...
            'optim_showcbar', true, ...
            'cmap_names', {{'discrim'}}, ...
            'cbar_locs', {{'east'}}, ...
            'demo_showyticklabel', true, ...
            'demo_xlim', [-58, -51], ...
            'optim_xlim', [0.4,1], ...
            'optim_dV_ylim', [-8,1], ...
            'optim_perf_ylim', [-0.05,1.05]), ...
    'recognition_TPR', ...
        struct(...
            'title', 'discrim$_{\mathrm{acc}}$ and recog$_{\mathrm{TPR}}$', ...
            'demo_showlegend', true, ...
            'optim_showlegend', false, ...
            'optim_showcbar', true, ...
            'cmap_names', {{'recognition_postradeoffs', 'recognition_TPR'}} , ...
            'cbar_locs', {{'west', 'east'}}, ...
            'demo_showyticklabel', false, ...
            'demo_xlim', [-63, -50], ...
            'optim_xlim', [0.3,1], ...
            'optim_dV_ylim', [-10,2], ...
            'optim_perf_ylim', [-0.05,1.05]) ...
);

recog_fields = fieldnames(recogplt_opts);


graphic_setdefault(15, ...
    'DefaultAxesMinorGridAlpha', 0.05, ...
    'DefaultAxesMinorGridLineStyle', '-', ...
    'DefaultTextInterpreter', 'latex', ...
    'DefaultLegendInterpreter', 'latex', ...
    'DefaultAxesLineWidth', 1.5, ...
    'DefaultLegendFontSize', 15, ...
    'DefaultLineLineWidth', 3);

figure('position', [0, 0, 0.7, 0.7]);
ax_list = tight_subplot(3, 2, ...
    [0.1, 0.05], ... gap
    [0.1,0.05], ... margin height
    [0.12, 0.02], ... margin width
    [0.3,1,0.6], ... ratio height
    [1,1], ... ratio width
    true, ... will plot 
    {'FontSize', 15, 'TitleFontSize', 1.3, 'LabelFontSize', 1.2}, ... axes options
    [1,0.2] ... gap height scale
);



for k = 1:length(recog_fields)
    field_recog = recog_fields{k};
    plt_opt = recogplt_opts.(field_recog);
    
    % demo plot
    set(gcf, 'CurrentAxes', ax_list(1,k));
    hold(gca, 'on');
    y_tick_lbls = cell(length(select_pairs), 1);
    
    for i = 1:length(select_pairs)
        select_alpha_val = select_pairs{i}.alpha;
        select_alpha_ind = find_nearest(alpha_W_vec, select_alpha_val, 'ind');
        select_alpha_val = alpha_W_vec(select_alpha_ind);
        
        select_incomp_val = select_pairs{i}.inp;
        select_incomp_ind = find_nearest(percent_complete_input_vec, select_incomp_val, 'ind');
        select_incomp_val = percent_complete_input_vec(select_incomp_ind);
        
        discrim_val = squeeze(Vm_optim.discrimination_accuracy(select_alpha_ind, select_incomp_ind, select_novlp_ind));
        recog_val = squeeze(Vm_optim.(field_recog)(select_alpha_ind, select_incomp_ind, select_novlp_ind));
        
        legend_on_tag = 'showlegend';
        if i > 1
            legend_on_tag = 'none';
        end
        
        plot([discrim_val, recog_val], [i, i], '-k', 'linewidth', 2)
        plot(discrim_val, i, ...
            'linestyle', 'none',  ...
            'marker', '>', ...
            'markersize', 12,  ...
            'markerfacecolor', cmaps.discrim(end,:), ...
            'markeredgecolor', 'none', ...
            'displayname', lgnd_shorts.discrim, ...
            'tag', legend_on_tag);
        
        plot(recog_val, i, ...
            'linestyle', 'none', ...
            'marker', '<', ...
            'markersize', 12, ...
            'markerfacecolor', cmaps.(field_recog)(end,:), ...
            'markeredgecolor', 'none', ...
            'displayname', lgnd_shorts.(field_recog), ...
            'tag', legend_on_tag);
        
        if plt_opt.demo_showyticklabel
            y_tick_lbls{i} = sprintf('(%.2f,%.1f)', select_alpha_val, select_incomp_val);
        end
        
    end
    
    xline(convert_obj.Vt_base, ':k', 'linewidth', 2, 'displayname', 'base', 'tag', 'showlegend');
    xlabel('$V_{T, \mathrm{optim}}$ (mV)');
    
    if plt_opt.demo_showyticklabel
        set(gca, ...
            'ytick', 1:length(select_pairs), ...
            'yticklabels', y_tick_lbls);
    else
        set(gca, 'ycolor', 'none');
    end
    
    despline;
    xlim(plt_opt.demo_xlim);
    ylim([0,length(select_pairs)+1]);
    
    if plt_opt.demo_showlegend
        lgnd = legend(...
            findobj(gca, 'tag', 'showlegend'), ...
            'Location', 'northwest', ...
            'NumColumns', 3, ...
            'box', 'on');
        lgnd.Position(2) = lgnd.Position(2) + 0.05;
    end
    
    % optim dV plot
    
    for j = 1:length(plt_optim_order)
        optim2plt = plt_optim_order{j};
        
        set(gcf, 'CurrentAxes', ax_list(j+1,k));
        hold(gca, 'on');
        for i = 1:length(ind_percent_inp)
            ind_sel = ind_percent_inp(i);
            inp_comp = percent_complete_input_vec(ind_sel)*100;
            ind_combos = {':',ind_sel,select_novlp_ind};
            switch lower(optim2plt)
                case 'dv'
                    discrim_datvec = smooth_fun(squeeze(dVm_optim.discrimination_accuracy(ind_combos{:})));
                    recog_datvec   = smooth_fun(squeeze(dVm_optim.(field_recog)(ind_combos{:})));
                case 'max'
                    discrim_datvec = squeeze(optim_results.discrimination_accuracy(:,ind_sel,select_novlp_ind));
                    recog_datvec   = squeeze(optim_results.(field_recog)(ind_combos{:}));
                otherwise
                    error('stop! figure it out!');
            end
            
            plot(alpha_W_vec, discrim_datvec, ...
                '-', 'linewidth', 3, ...
                'color', cmaps.discrim(i,:), ...
                'tag', 'showlegend', ...
                'displayname', sprintf('discrim (%02d\\%%)', inp_comp));
            plot(alpha_W_vec, recog_datvec, ...
                '-','linewidth', 3, ...
                'color',  cmaps.(field_recog)(i,:), ...
                'tag', 'showlegend', ...
                'displayname', sprintf('recog (%02d\\%%)', inp_comp));
            
        end
        
        if j == length(plt_optim_order)
            xlabel('$\alpha_W$');
        else 
            set(gca, 'xcolor', 'none');
        end 
        
        xlim(plt_opt.optim_xlim);
        
        despline;
        switch lower(optim2plt)
            case 'dv'
                
                title(plt_opt.title);
                if k == 1
                    ylabel('$\Delta V_{T, \mathrm{optim}}$ (mV)');
                end
                ylim(plt_opt.optim_dV_ylim);
            case 'max'
                if k == 1
                    ylabel('max performance')
                end
                ylim(plt_opt.optim_perf_ylim);
        end
        
            alpha_colors = flipud(gray(length(select_pairs)+1)) * 0.85;
            alpha_colors = alpha_colors(2:end,:);
            for i = 1:length(select_pairs)
                select_alpha_val = select_pairs{i}.alpha;
                select_alpha_ind = find_nearest(alpha_W_vec, select_alpha_val, 'ind');
                select_alpha_val = alpha_W_vec(select_alpha_ind);
                xline(select_alpha_val, '--','linewidth', 2, 'color', [alpha_colors(i,:), 0.8]);
            end
        
        
        if j == 2 && plt_opt.optim_showlegend
            legend(...
                findobj(gca, 'tag', 'showlegend'), ...
                'location', 'southeast', ...
                'NumColumns', 2, ...
                'color', [1,1,1]*0.97, ...
                'box', 'on');
        end
        
        if j == 1 && plt_opt.optim_showcbar
            cmap_names = plt_opt.cmap_names; 
            cbar_locs = plt_opt.cbar_locs; 
            for i_cbar = 1:length(cmap_names)
                cmap_name = cmap_names{i_cbar};
                cbar_loc = cbar_locs{i_cbar};
                
                if i_cbar == 1
                    ax_cur = gca;
                else
                    ax_cur = axes(gcf, 'position', get(gca, 'position'), ...
                        'units', 'normalized', 'visible', 'off');
                end
                
                colormap(ax_cur, cmaps.(cmap_name)); 
                cbar = colorbar(ax_cur,...
                    'Location', cbar_loc, ...
                    'FontSize', 15, ...
                    'TickLength', 0, ...
                    'Box', 'on', ...
                    'Linewidth', 1);
                caxis(ax_cur,[0,1]);
                
                cbar.Position = cbar.Position .* [1,1,0.8,0.4] + [0,0,0,0];
                
                if k == 1
                    xlabel(cbar, '\% complete', 'FontSize', 18, 'Interpreter', 'latex');
                end
                
                cbar_labels = arrayfun(@(x) sprintf('%.1f', x), percent_complete_input_vec(ind_percent_inp), 'uni', 0);
                cbar_ticks = linspace(0,1,length(cbar_labels)+1);
                cbar_ticks = 0.5 * (cbar_ticks(1:end-1) + cbar_ticks(2:end));
                
                title(cbar, ['\textit{' lgnd_shorts.(cmap_name) '}'], 'FontSize', 15, 'Interpreter', 'latex');
                
                set(cbar, ...
                    'ticks', cbar_ticks, ...
                    'ticklabels', cbar_labels)
            end

            
        end
    end
end


fig_filename = fullfile(fig_path, 'Fig6b');
export_fig(fig_filename, '-pdf', '-p0.02');
close; 