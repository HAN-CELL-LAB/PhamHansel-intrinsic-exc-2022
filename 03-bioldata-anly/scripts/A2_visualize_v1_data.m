% THIS FILE WILL NOT RUN WITHOUT V1 DATA 
% SEE README FILE 

clc; clear; close all; 
run start_up.m; 

data_path = 'data/v1-data/v1-invivo_Li2020.mat';
fig_path = fullfile('figures/v1-data');
if ~exist(fig_path, 'dir')
    mkdir(fig_path);
end

load(data_path);

%% Visualize violin
mean_Vrest =  mean(data.Vrest.data);

combine_name_and_unit = @(s) sprintf('$%s (%s)$', s.latex, s.unit);

violin_color = [0.4,0.4,0.4];
scatterviolin_color = [0,0,0,0.2];

common_violinopts = {...
    'ShowMean', true,  ....
    'ShowNotches', false, ...
    'ViolinColor', violin_color};

fig_description = 'V1 in vivo (Li et al 2020)';

figure('DefaultAxesTickLabelInterpreter', 'latex',  ...
    'units', 'normalized', ...
    'Position', [0,0,0.5,0.95], ...
    'DefaultLineLineWidth', 2, ...
    'DefaultAxesFontSize', 30, ...
    'DefaultAxesLabelFontSize', 1.2, ...
    'DefaultAxesLineWidth', 2);

subplot(1,3,[1,2]); hold on;
MP_prop_vars = {'Vrest', 'Vthres'};
MP_prop_data = cell2mat(cellfun(@(x)  data.(x).data, MP_prop_vars, 'uni', 0));
MP_prop_name = {'rest ($V_R$)', 'thres ($V_T$)'};

vp = violinplot(MP_prop_data, MP_prop_name, common_violinopts{:});

vp_scatters = [vp.ScatterPlot];
arrayfun(@(x) set(x, 'SizeData', 150),  vp_scatters);

plot(vertcat(vp_scatters.XData), vertcat(vp_scatters.YData), 'color', scatterviolin_color, 'linewidth', 1);
xlim([0.5,2.5]);

ylabel('membrane potential (mV)');
despline;
horz_pos_ax(0.02, 0.98);

subplot(1,3,3); hold on;
vp = violinplot(data.dV.data, {combine_name_and_unit(data.dV)}, common_violinopts{:});

vp_scatters = [vp.ScatterPlot];
arrayfun(@(x) set(x, 'SizeData', 150),  vp_scatters);
xlim([0.5,1.5]);
ylim([0, 40]);
xlabel('$\Delta V_{TR} = V_T - V_R$');
ylabel('(mV)'); 
hide_only_axis('x');
horz_pos_ax(0.02);


pause(1);
fig_name = fullfile(fig_path, 'v1-violin.pdf');
exportgraphics(gcf, fig_name); 
close;

%% Helper
function horz_pos_ax(dx, kw)
if ~exist('kw', 'var'), kw = 1; end 
pos = get(gca, 'position');
pos(1) = pos(1) + dx; 
pos(3) = pos(3) * kw; 
set(gca, 'position', pos)
end