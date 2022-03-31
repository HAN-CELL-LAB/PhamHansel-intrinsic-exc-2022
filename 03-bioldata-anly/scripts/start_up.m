%% Example start up scripts to call before every script

% Add path
% addpath(genpath('.')); 
addpath(genpath('functions'));
addpath(genpath('extpkg'));

% Remove path
% rmpath(genpath('./scratch'));
% rmpath(genpath('./figures'));
% rmpath(genpath('./docs'));
% rmpath(genpath('./demo'));

% Reset default graphic options 
graphic_setdefault(15, ...
    'DefaultAxesMinorGridAlpha', 0.05, ...
    'DefaultAxesMinorGridLineStyle', '-', ...
    'DefaultTextInterpreter', 'latex', ...
    'DefaultLegendInterpreter', 'latex', ...
    'DefaultStemMarkerSize', 1, ...
    'DefaultStemlineWidth', 1.5);

% Any additional settings
% set(0, 'DefaultFigureWindowStyle','docked');
