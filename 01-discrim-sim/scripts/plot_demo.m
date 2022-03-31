clc; clear; close all; 

addpath(genpath('extpkg'));
addpath(genpath('functions'));
graphic_setdefault(15, ...
    'DefaultAxesMinorGridAlpha', 0.05, ...
    'DefaultAxesMinorGridLineStyle', '-', ...
    'DefaultTextInterpreter', 'latex', ...
    'DefaultLegendInterpreter', 'latex', ...
    'DefaultStemMarkerSize', 1, ...
    'DefaultStemlineWidth', 1.5, ...
    'DefaultFigureWindowStyle','normal');

%%
Jd = [2/8,3/8,4/7,5/7,4/5];
Hd = [2/8,3/8,4/8,5/8,4/8];

Jd_color = return_colorbrewer('Reds',5);
Jd_color = Jd_color(end,:);

Hd_color = return_colorbrewer('Blues',5);
Hd_color = Hd_color(end,:);

common_styles = {'linewidth', 2.5, 'marker', '.', 'markersize', 40};
figure('units','normalized','Position',[0.1,0.1,0.1,0.6]); hold on;
plot(Jd, 'color', Jd_color, common_styles{:});
plot(Hd, 'color', Hd_color, common_styles{:});
view(90,90);

set(gca,'xcolor','none','ycolor','none');

exportgraphics(gcf, 'misc/demo.pdf', 'ContentType', 'vector')