clc; clear; close all; 
run start_up.m; 

fig_path = fullfile('figures/conversion');
if ~exist(fig_path, 'dir')
    mkdir(fig_path);
end

load('data/V1invivo_Li2020.mat');


%% Visualize

mean_Vrest =  mean(data.Vrest.data);

combine_name_and_unit = @(s) sprintf('$%s (%s)$', s.latex, s.unit);

violin_color = [0.4,0.4,0.4];
scatterviolin_color = [0,0,0,0.2];


common_violinopts = {'ShowMean', true, 'ShowNotches', true, 'ViolinColor', violin_color};

fig_description = 'V1 in vivo (Li et al 2020)';

figure('DefaultAxesTickLabelInterpreter', 'latex',  ...
    'units', 'normalized', 'Position', [0,0,1,1], ...
    'DefaultLineLineWidth', 2, ...
    'DefaultAxesFontSize', 30, ...
    'DefaultAxesLabelFontSize', 1.2, ...
    'DefaultAxesLineWidth', 1.5);

annotation('textbox', 'String', fig_description, ...
    'FontSize', 30, 'Interpreter', 'latex', 'LineStyle', 'none', ...
    'Position', [0, 0.92, 0.9, 0.08], ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');

subplot(2,3,[1,4]); hold on;
MP_prop_vars = {'Vrest', 'Vthres'};
MP_prop_data = cell2mat(cellfun(@(x)  data.(x).data, MP_prop_vars, 'uni', 0));
MP_prop_name = {'rest ($V_R$)', 'thres ($V_T$)'};

vp = violinplot(MP_prop_data, MP_prop_name, common_violinopts{:});

vp_scatters = [vp.ScatterPlot];
arrayfun(@(x) set(x, 'SizeData', 50),  vp_scatters);

plot(vertcat(vp_scatters.XData), vertcat(vp_scatters.YData), 'color', scatterviolin_color, 'linewidth', 1);
% set(get(gca, 'XAxis'), 'TickLength', [0,0]);
xlim([0.5,2.5]);

ylabel('membrane potential (mV)');
despline;

subplot(2,6,[3,9]); hold on;
vp = violinplot(data.dV.data, {combine_name_and_unit(data.dV)}, common_violinopts{:});

vp_scatters = [vp.ScatterPlot];
arrayfun(@(x) set(x, 'SizeData', 50),  vp_scatters);
xlim([0.5,1.5]);

xlabel('$\Delta V_{TR} = V_T - V_R$');
ylabel('(mV)'); 
hide_only_axis('x');

ax1 = subplot(2,3,3); hold on;
plot(data.Vrest.data, data.Vthres.data, '.k', 'markersize', 30, 'tag', 'sharexaxis');
xline(mean(data.Vrest.data), ':k', 'linewidth', 2);
id_vals = [-70, -45];
plot(id_vals, id_vals, '--k', 'linewidth', 2);

xlabel('$V_R$ (mV)'); ylabel('$V_T$ (mV)'); 
xlim([-80, -40]); despline; 

ax2 = subplot(2,3,6); hold on;
plot(data.Vrest.data, data.dV.data, '.k', 'markersize', 30,  'tag', 'sharexaxis');
xlabel('$V_R$ (mV)'); ylabel('$\Delta V_{TR}$ (mV)'); 
xlim([-80, -40]); despline; 

ax1.Position(1) = ax1.Position(1) - 0.07;
ax2.Position(1) = ax1.Position(1);

fig_name = fullfile(fig_path, 'Li2020-data');
export_fig(fig_name,  '-r300', '-p0.02');
pause(0.5); close;
