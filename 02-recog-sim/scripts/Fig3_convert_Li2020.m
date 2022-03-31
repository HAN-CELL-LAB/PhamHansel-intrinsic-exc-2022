clc; clear; close all;
run start_up;

fig_path = fullfile('figures/conversion');
if ~exist(fig_path, 'dir')
    mkdir(fig_path);
end

load('data/V1invivo_Li2020.mat');

%% Define common variables
theta_range = [0, 1];
theta_base = 0.5;
num_points = 10001;

pop_mean_Vr = mean(data.Vrest.data);
pop_mean_Vt = mean(data.Vthres.data);

%% Convert using whole range of whole population 
convert_data_obj = struct;

source_Vt_range = [min(data.Vrest.data), max(data.Vthres.data)];
source_Vt_base = mean(data.Vthres.data);
get_Vt_max = @(vr, vt) vr + diff(source_Vt_range)*(vt-vr)/(source_Vt_base-source_Vt_range(1));

convert_data_obj = struct(...
    'description', 'individual, thres as base', ...
    'latex',  ['$[V_R^{(i)}, V_T^{(i)}] \rightarrow [\theta_{min}^{(i)}, \theta_{base}^{(i)}]$' ...
    newline '$V_{T,max}^{(i)} = $ extrapolation' ], ...
    'Vt_range', [data.Vrest.data, get_Vt_max(data.Vrest.data, data.Vthres.data)], ...
    'Vt_base', data.Vthres.data);

fontsize_latex = 10; 

modify_latex_fontsize = @(s) regexprep(s, '^',...
    sprintf('\\\\fontsize{%f}{0}\\\\selectfont', fontsize_latex), ...
    'emptymatch', 'lineanchors'); 
text_align_latex = @(s) regexprep(s, {'^', '$', newline}, ...
    {'\\begin{tabular}{l} ', ' \\end{tabular}', '\\\\'}, ...
    'emptymatch'); 
modify_style_latex = @(s) text_align_latex(modify_latex_fontsize(s)); 

convert_data_obj = setfield(convert_data_obj, 'latex', modify_style_latex(convert_data_obj.latex));
convert_wholepop_funlist = {@convert_Vt_to_theta_quadratic};

%% Plotting 
figure('DefaultAxesTickLabelInterpreter', 'latex',  ...
    'units', 'normalized', 'Position', [0,0,0.8,0.7], ...
    'DefaultLineLineWidth', 2, 'DefaultAxesFontSize', 30,...
    'DefaultAxesLabelFontSize', 1,  'DefaultAxesTitleFontSize', 1);


ax1 = subplot(1,2,1); hold(ax1, 'on');
ax2 = subplot(1,2,2); hold(ax2, 'on');

convert_fun = @convert_Vt_to_theta_quadratic;

indiv_convert_results = arrayfun(@(x) ...
    convert_fun(convert_data_obj.Vt_range(x,:), convert_data_obj.Vt_base(x), theta_range, theta_base, num_points), ...
    1:size(convert_data_obj.Vt_range,1), 'uni', 1);

tmp_res = indiv_convert_results(1);

fieldpairs_to_plot = {{'Vt_vec', 'theta_vec'}, {'delta_Vt_vs_base', 'delta_theta_vs_base'}};
line_plt_opts_eachcell = {'Color', [indiv_convert_results(1).color, 0.1], 'LineWidth', 2};
line_plt_opts_meanpop = {'Color', [indiv_convert_results(1).color, 0.8], 'LineWidth', 3, 'DisplayName', 'mean'};
axes_for_plot = [ax1, ax2];
n_sample = length(indiv_convert_results);

for j = 1:length(fieldpairs_to_plot)
    fieldpair = fieldpairs_to_plot{j};
    arrayfun(@(res) plot(axes_for_plot(j), res.(fieldpair{1}), res.(fieldpair{2}), line_plt_opts_eachcell{:}), indiv_convert_results);
    
    mean_of_fields = cellfun(@(s) mean(vertcat(indiv_convert_results.(s))), fieldpair, 'uni', 0);
    errbar_of_fields = cellfun(@(s) std(vertcat(indiv_convert_results.(s))), fieldpair, 'uni', 0);
    
    plot(axes_for_plot(j), mean_of_fields{1}, mean_of_fields{2}, line_plt_opts_meanpop{:}, 'tag', 'legendon');
    
    fill(axes_for_plot(j), ...
        [mean_of_fields{1}+errbar_of_fields{1}, fliplr(mean_of_fields{1}-errbar_of_fields{1})], ...
        [mean_of_fields{2}, fliplr(mean_of_fields{2})], ...
        indiv_convert_results(1).color, 'LineStyle', 'none', 'FaceAlpha', 0.1, ...
        'DisplayName', 'SD($x$)', 'Tag', 'legendon');
    
    tmp_res.(fieldpair{1}) = mean_of_fields{1};
    tmp_res.(fieldpair{2}) = mean_of_fields{2};
    
end

tmp_res.Vt_rangeNbase = mean(vertcat(indiv_convert_results.Vt_rangeNbase));
tmp_res.Vt_range = tmp_res.Vt_rangeNbase([1,end]);
tmp_res.Vt_base = mean(vertcat(indiv_convert_results.Vt_base));

convert_data_obj.results.(tmp_res.description) = tmp_res;

xline(ax1, pop_mean_Vr, '--k', 'linewidth', 2, ...
    'displayname', '$\overline{V_R}$', 'tag', 'legendon');
xline(ax1, pop_mean_Vt, '--', 'color', 0.7*[1,1,1], 'linewidth', 2, ...
    'displayname', '$\overline{V_T}$', 'tag', 'legendon');
yline(ax1, theta_base, ':k', 'linewidth', 1);
xline(ax2, 0, ':k', 'linewidth', 1);
yline(ax2, 0, ':k', 'linewidth', 1);


title(ax1, convert_data_obj.latex);
xlabel(ax1, '$x \sim V_T$ (mV)');
ylabel(ax1, '$y \sim \theta$');

xlabel(ax2, '$V_T - V_{T,b}$ (mV)');
ylabel(ax2, '$\theta - \theta_b$');
title(ax2, 'converted change');

despline([ax1, ax2]);
ylim(ax1, [-0.03,1.03]);
ylim(ax2, [-0.03,1.03]-0.5);
xlim(ax2, [-25, 20]);
lgnd = legend(ax1, findobj(ax1, 'tag', 'legendon'), 'fontsize', 20);
title(lgnd, tmp_res.latex);


% fig_name = fullfile(fig_path, 'conversion-Li2020-individualcurved');
% export_fig(fig_name,  '-r200', '-p0.02');
% pause(0.5); close;

% 
% save('data/conversion_Li2020.mat', ...
%     'convert_data_obj', 'theta_range', 'theta_base', 'num_points');


%% Functions for conversion 
function res = convert_Vt_to_theta_linear(Vt_range, ~, theta_range, theta_base, num_points)

p1 = polyfit(Vt_range, theta_range, 1);
Vt_vec = linspace(Vt_range(1), Vt_range(2), num_points);
theta_vec = polyval(p1, Vt_vec);
Vt_base = Vt_vec(find_nearest(theta_vec, theta_base, 'ind'));

res = struct;
res.color = [0,0,0];
res.description = 'linear';
res.latex = '$y \sim ax + b$';
res.Vt_vec = Vt_vec;
res.theta_vec = theta_vec;
res.Vt_base = Vt_base;
res.delta_Vt_vs_base = Vt_vec - Vt_base;
res.delta_theta_vs_base = theta_vec - theta_base;

res.Vt_rangeNbase = sort([Vt_range, Vt_base]);
res.theta_rangeNbase = sort([theta_range, theta_base]);

end

function res = convert_Vt_to_theta_quadratic(Vt_range, Vt_base, theta_range, theta_base, num_points)
res = struct;
% res.color = [0,0,1];
res.color = [0.1,0.4,0.2];

res.description = 'quadratic';
res.latex = '$y \sim ax^2 + bx + c$';

res.Vt_rangeNbase = sort([Vt_range, Vt_base]);
res.theta_rangeNbase = sort([theta_range, theta_base]);

p2 = polyfit(res.Vt_rangeNbase, res.theta_rangeNbase, 2);
Vt_vec = linspace(Vt_range(1), Vt_range(2), num_points);
theta_vec = polyval(p2, Vt_vec);

res.fit_fun = @(x) polyval(p2, x); 
res.Vt_vec = Vt_vec;
res.theta_vec = theta_vec;

res.delta_Vt_vs_base = Vt_vec - Vt_base;
res.delta_theta_vs_base = theta_vec - theta_base;

res.Vt_base = Vt_base;
end


function res = convert_Vt_to_theta_flipquad(Vt_range, Vt_base, theta_range, theta_base, num_points)
res = struct;
res.color = [1,0,0];
res.description = 'flipped_quadratic';
res.latex = '$x \sim ay^2 + by + c$';

res.Vt_rangeNbase = sort([Vt_range, Vt_base]);
res.theta_rangeNbase = sort([theta_range, theta_base]);

p2 = polyfit(res.theta_rangeNbase, res.Vt_rangeNbase, 2);
theta_vec = linspace(theta_range(1), theta_range(2), num_points);
Vt_vec = polyval(p2, theta_vec);

res.fit_fun = @(x) polyval(p2, x); 
res.Vt_vec = Vt_vec;
res.theta_vec = theta_vec;

res.delta_Vt_vs_base = Vt_vec - Vt_base;
res.delta_theta_vs_base = theta_vec - theta_base;

res.Vt_base = Vt_base;
end

