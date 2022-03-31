t = 0:0.01:30;
tau_1 = 0.5; 
tau_2 = 7; 

syn = -exp(-t/tau_1) + exp(-t/tau_2);
syn = syn/max(syn);

figure; hold on; 
plot(t, syn, '-k', 'linewidth', 10);
plot(t, 2*syn, '-b', 'linewidth', 10);
plot(t, 0.5*syn, '-r', 'linewidth', 10);

text(1, 1.2, 'LTP', 'Color', 'b', 'fontsize', 70);
text(1, 0.2, 'LTD', 'Color', 'r', 'fontsize', 70);

set(gca, 'xcolor', 'none', 'ycolor', 'none'); 
title('EPSP', 'fontsize', 80);
% export_fig('figures/demo_sp', '-r200', '-p0.02');

%%
x = -8:0.01:8; 

figure; hold on; 
plot(x, 1./(1+exp(-x)), '-k', 'linewidth', 10);
plot(x, 1./(1+exp(-(x - 2))), '-r', 'linewidth', 10);
plot(x, 1./(1+exp(-(x + 2))), '-b', 'linewidth', 10);

text(-4, 0.8, '$\uparrow$ IE', 'Color', 'b', 'fontsize', 70, 'Interpreter', 'latex');
text(2, 0.2, '$\downarrow$ IE', 'Color', 'r', 'fontsize', 70, 'Interpreter', 'latex');
set(gca, 'xcolor', 'none', 'ycolor', 'none'); 
title('f-I curve', 'fontsize', 80); 

% export_fig('figures/demo_ip', '-r200', '-p0.02');

%%
x = -10:0.01:10; 
g = 1000; 
figure; hold on; 
plot(x, 1./(1+exp(-g*x)), '-k', 'linewidth', 10);
plot(x, 1./(1+exp(-g*(x - 5))), '-r', 'linewidth', 10);
plot(x, 1./(1+exp(-g*(x + 5))), '-b', 'linewidth', 10);

text(-3.5, 0.85, '$\uparrow$ IE', 'Color', 'b', 'fontsize', 80, 'Interpreter', 'latex');
text(-3.5, 0.65, '$\downarrow \theta$', 'Color', [0.8,0.8,1], 'fontsize', 80, 'Interpreter', 'latex');

text(1, 0.2, '$\downarrow$ IE', 'Color', 'r', 'fontsize', 80, 'Interpreter', 'latex');
text(1, 0.4, '$\uparrow \theta$', 'Color', [1,0.8,0.8], 'fontsize', 80, 'Interpreter', 'latex');

annotation('textarrow',[0.5, 1/3],[0.5, 0.5], ...
    'Color', [0.8,0.8,1], 'LineWidth', 5, ...
    'HeadLength', 50, 'HeadWidth', 50);

annotation('textarrow',[0.535, 0.7],[0.5, 0.5], ...
    'Color', [1,0.8,0.8], 'LineWidth', 5, ...
    'HeadLength', 50, 'HeadWidth', 50);


set(gca, 'xcolor', 'none', 'ycolor', 'none'); 
title('activation', 'fontsize', 80); 

% export_fig('figures/demo_ip_heaviside', '-r200', '-p0.02');

%%
ip_color = [0.7,0.7,0.7,0.5]; 

x = -10:0.01:10; 
g = 1000; 
figure; hold on; 
plot(x, 1./(1+exp(-g*x)), '-k', 'linewidth', 10);
plot(x, 1./(1+exp(-g*(x - 5))), '-', 'color', ip_color,  'linewidth', 10);
plot(x, 1./(1+exp(-g*(x + 5))), '-', 'color', ip_color, 'linewidth', 10);

text(-3.5, 0.85, '$\uparrow$ IE', 'Color', ip_color, 'fontsize', 80, 'Interpreter', 'latex');
text(-3.5, 0.65, '$\downarrow \theta$', 'Color', ip_color, 'fontsize', 80, 'Interpreter', 'latex');

text(1, 0.2, '$\downarrow$ IE', 'Color', ip_color, 'fontsize', 80, 'Interpreter', 'latex');
text(1, 0.4, '$\uparrow \theta$', 'Color', ip_color, 'fontsize', 80, 'Interpreter', 'latex');

annotation('textarrow',[0.5, 1/3],[0.5, 0.5], ...
    'Color', ip_color, 'LineWidth', 5, ...
    'HeadLength', 50, 'HeadWidth', 50);

annotation('textarrow',[0.535, 0.7],[0.5, 0.5], ...
    'Color', ip_color, 'LineWidth', 5, ...
    'HeadLength', 50, 'HeadWidth', 50);


set(gca, 'xcolor', 'none', 'ycolor', 'none'); 
title('activation', 'fontsize', 100); 

export_fig('figures/demo_ip_heaviside', '-r200', '-p0.02');