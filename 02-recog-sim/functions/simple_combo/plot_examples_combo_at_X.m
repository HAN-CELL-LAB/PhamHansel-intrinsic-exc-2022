function res = plot_examples_combo_at_X(N, num_thres, example_id, distrib_info, willplot, fig_subpath, show_figdescription, exportfig_type)
if ~exist('willplot', 'var'), willplot = false; end 
if ~exist('show_figdescription', 'var'), show_figdescription = true; end 
if ~exist('fig_subpath', 'var'), fig_subpath = ''; end 
if ~exist('exportfig_type', 'var'), exportfig_type = '-r300'; end

distrib_type = distrib_info.type;
switch lower(distrib_type)
    case 'lognorm'
        sigma_lognorm = distrib_info.param;
        genWfun = @(varargin) lognrnd(0,sigma_lognorm,varargin{:});
        distrib_str = sprintf('$\\sigma_{\\mathrm{lognorm}} = %.2f$', sigma_lognorm);
        fig_name = sprintf('combo-at-X_lognorm-sigma-%.2f_example%02d', sigma_lognorm, example_id); 
    case 'equal'
        genWfun = @(varargin) ones(varargin{:});
        distrib_str = 'equal weights';
        fig_name = sprintf('combo-at-X_equal_example%02d', example_id); 
    case 'uniform'
        genWfun = @(varargin) rand(varargin{:});
        distrib_str = 'weights from uniform distribution before normalized';
        fig_name = sprintf('combo-at-X_uniform_example%02d', example_id); 
    otherwise
        error('"%s" as type is not allowed', distrib_type);
end
description_str = '$Y=1$ if $\sum_i W_i x_i > \theta$';

fig_description =  sprintf(['N = %d ' newline '%s' newline '%s'], N, description_str, distrib_str);

W = genWfun(1,N);
W = W/sum(W);
W = sort(W, 'descend');

thres = linspace(0,1,num_thres);
X = de2bi(0:(2^N-1));

Ypre = W * X';
[Ypre, ind_sort] = sort(Ypre, 'descend');
X = X(ind_sort, :);

Yact = Ypre' > thres;
nYact = sum(Yact, 1);

res = struct();
res.W = W; 
res.nYact = nYact;
res.Ypre = Ypre; 
res = structfun(@to_row_vec, res, 'uni', 0);

if ~willplot
    return; 
end

figure('units','normalized','position',[0,0.0,0.65,1]);
if show_figdescription
    annotation(gcf,'textbox',...
        'units','normalized', 'position',[0.02, 0.8, 0.3, 0.15], ...
        'FitBoxToText', 'off', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'cap', ...
        'LineStyle', 'none', 'FontSize', 20, 'Interpreter', 'latex', ...
        'String', fig_description);
end
ax_weight = axes('units', 'normalized', 'position', [0.26    0.80    0.25    0.15]);
ax_thres  = axes('units', 'normalized', 'position', [0.62    0.80    0.25    0.15]);
ax_combo  = axes('units', 'normalized', 'position', [0.26    0.10    0.25    0.64]);
ax_spike  = axes('units', 'normalized', 'position', [0.62    0.10    0.25    0.64]);
ax_sumW   = axes('units', 'normalized', 'position', [0.10    0.10    0.11    0.64]);

linkaxes([ax_sumW, ax_combo, ax_spike], 'y');
linkaxes([ax_combo, ax_weight], 'x');

hold(ax_combo, 'on'); colormap('gray');
image_with_strict_limits(ax_combo, X');

hold(ax_spike, 'on'); colormap('gray');
image_with_strict_limits(ax_spike, Yact');

hold(ax_weight, 'on');
axwb = bar(ax_weight, W,  'EdgeAlpha', 0.8, 'FaceAlpha', 0.1, 'FaceColor', 'k', 'LineWidth', 1.2);
axwb.BaseLine.LineStyle = 'none';
ylim(ax_weight,[0,0.3]);

hold(ax_thres, 'on');
axth = bar(ax_thres, thres, nYact, 1,  'EdgeAlpha', 0, 'FaceAlpha', 0.1, 'FaceColor', 'k');
axth.BaseLine.LineStyle = 'none';
hold(ax_sumW, 'on');
barh(ax_sumW, Ypre, 1, 'EdgeAlpha', 0, 'FaceAlpha', 0.1, 'FaceColor', 'k')

set(ax_sumW, 'XDir', 'reverse');
ylabel(ax_sumW, {'sorted weighted sum of inputs',...
    '(combinations sorted by $\sum W_i x_i$)'});
hide_only_axis(ax_sumW,'xy');

set(ax_combo, 'box', 'on', 'xtick', '', 'ytick', '');
xlabel(ax_combo, {'combination of inputs', '(white = on)'});

set(ax_spike, 'box', 'on', 'xtick', '', 'ytick', '');

title(ax_weight, 'sorted inputs weights');
xticks(ax_weight, 1:N);
set(ax_weight, 'xcolor', 'none');

xlabel(ax_spike, {'threshold $\theta$','(white $\sim Y = 1$)'});
title(ax_thres, '\# of comb. with $Y=1$');

export_fig(fullfile('figures', fig_subpath, fig_name), exportfig_type, '-p0.02'); 
close; 
end