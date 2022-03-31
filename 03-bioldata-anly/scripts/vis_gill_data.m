clc; clear; close all;
run start_up.m

graphic_setdefault(15, ...
    'DefaultAxesMinorGridAlpha', 0.05, ...
    'DefaultAxesMinorGridLineStyle', '-', ...
    'DefaultTextInterpreter', 'latex', ...
    'DefaultLegendInterpreter', 'latex', ...
    'DefaultAxesTitleFontSize', 1.1, ...
    'DefaultAxesLabelFontSize', 1.1, ...
    'DefaultAxesLineWidth', 2, ...
    'DefaultLineLineWidth', 3);

%% Define paths
main_data_path = 'data/gill-data/Layer II.III S1 BC IP';

fig_path = 'figures/gill-data/Layer II.III S1 BC IP'; 
if ~exist(fig_path, 'dir')
    mkdir(fig_path); 
end

INFO_path = fullfile(main_data_path, 'info');
PROCESS_path = fullfile(main_data_path, 'processed');

%% Load data
load(fullfile(PROCESS_path, 'analysis-summary.mat'));

analysis_results = analysis_table.analysis;
representative_traces = analysis_table.representative;

cell_fullIDs = arrayfun(@(x) x.full_ID, analysis_table.info, 'uni', 0); 
exp_groups = arrayfun(@(x) x.Group, analysis_table.info, 'uni', 0);

%% Defining IDs and figure names 
cellids_latex = cellfun(@(x,y) sprintf('\\textbf{%s} (\\textit{%s})',x,y), cell_fullIDs, exp_groups, 'uni', 0);
fignames = cellfun(@(x,y) sprintf('CELL-%s--GROUP-%s',x,y), cell_fullIDs, exp_groups, 'uni', 0);

%% Selection
condition_source = exp_groups;
select_conditions = contains(condition_source, {'WT'});

select_analysis = analysis_results(select_conditions);
select_representative = representative_traces(select_conditions);
select_cellidlatex = cellids_latex(select_conditions); 
select_fignames = fignames(select_conditions); 

select_num_cells = length(select_analysis);

%% Define generals 

color_1 = [0.8,0.3,0.3,0.8];
color_2 = [0.3,0.3,0.9,0.8];
color_3 = [0.25,0.5,0.2,0.8];

repr_colors = [0.7,0.7,0.7;0.7,0.7,0.95;0.2,0.2,1];
repr_lgdn = {'base', 'post 20\%', 'post 90\%'};
repr_lbls_gen = struct('x', 'time (ms)', 'y', '$V_m$ (mV)', 'ttl', {{'representative (spikes)', 'representative ($R_{in}$)'}});

t0 = tic; 

for i = 1:select_num_cells
    cell_id = select_cellidlatex{i}; 
    anly_struct = select_analysis(i);
    repr_struct = select_representative(i);
    stim_info = anly_struct.src(1).configs.stim_info;
    anly_xvec = anly_struct.time_vec;
    
    repr_lbls = repr_lbls_gen; 
    repr_lbls.ttl{1} = sprintf('%s - %s',  cell_id, repr_lbls_gen.ttl{1});
    figure;
    ax1 = subplot(3,3,[1,2]);
    ax2 = subplot(3,3,3);
    plot_representative_traces([ax1,ax2],repr_struct,stim_info,repr_colors,repr_lgdn,repr_lbls)
    despline([ax1,ax2]);
    
    subplot(334); hold on; 
    yyaxis left; 
    plot_with_smoothing(anly_xvec,anly_struct.num_spikes,...
        struct('x','trace \#','y','\# spikes','ttl','baseline (trace $< 0$) vs post (trace $> 0$)'), color_1);
    yyaxis right; 
    plot_with_smoothing(anly_xvec,anly_struct.latency,...
        struct('x','trace \#','y','latency (ms)','ttl',''), color_2);
 
    subplot(335); hold on; 
    yyaxis left; 
    plot_with_smoothing(anly_xvec,anly_struct.Rin,...
        struct('x','trace \#','y','$R_{in}$','ttl',''),color_1);
    
    yyaxis right; 
    plot_with_smoothing(anly_xvec,anly_struct.Vthres_first - anly_struct.Vrest_1,...
        struct('x','trace \#','y','$V_{\mathrm{thres}}^{(1)} - V_{\mathrm{rest-1}}$','ttl',''),color_2);
    
    subplot(336); hold on;
    yyaxis left;
    plot_with_smoothing(anly_xvec,anly_struct.AP_Vm_mean - anly_struct.Vbase,...
        struct('x','trace \#','y','$V_{\mathrm{AP,mean}} - V_{\mathrm{pre-pulse}}$','ttl',''),color_1);
    
    yyaxis right;
    plot_with_smoothing(anly_xvec,anly_struct.fAHP_Vm_mean - anly_struct.Vbase,...
        struct('x','trace \#','y','$V_{\mathrm{fAHP,mean}} - V_{\mathrm{pre-pulse}}$','ttl',''), color_2);

    subplot(6,3,13); hold on;
    
    plot_with_smoothing(anly_xvec,anly_struct.Vthres_first,...
        struct('x','','y','(mV)','ttl',''),color_2,true,'$V_{\mathrm{thres}}^{(1)}$ (first thres)');
    
    plot_with_smoothing(anly_xvec,anly_struct.Vrest_1,...
        struct('x','','y','(mV)','ttl',''),color_1,true,'$V_{\mathrm{rest-1}}$ (Tuan)');
    
    set(gca,'xcolor', 'none','ycolor','k'); 
    legend_on_selected(gca);

    subplot(6,3,14); hold on;
   
    plot_with_smoothing(anly_xvec,anly_struct.Vthres_mean_afterfirst,...
        struct('x','','y','(mV)','ttl',''),color_2,true,'$V_{\mathrm{thres,mean}}^{(> 1)}$ (mean after first)');
    
    plot_with_smoothing(anly_xvec,anly_struct.Vrest_2,...
        struct('x','','y','(mV)','ttl',''),color_1,true,'$V_{\mathrm{rest-2}}$ (Christian)');
    set(gca,'xcolor', 'none','ycolor','k'); 
    
    legend_on_selected(gca);

    subplot(6,3,15); hold on;
   
    plot_with_smoothing(anly_xvec,anly_struct.Vthres_mean_all,...
        struct('x','','y','(mV)','ttl',''),color_2,true,'$V_{\mathrm{thres,mean}}^{(all)}$ (mean all)');
    
    plot_with_smoothing(anly_xvec,anly_struct.Vbase,...
        struct('x','','y','(mV)','ttl',''),color_1,true,'$V_{\mathrm{pre-pulse}}$');
    set(gca,'xcolor', 'none','ycolor','k'); 
    
    legend_on_selected(gca);
    
    subplot(6,3,16); hold on;
   
    plot_with_smoothing(anly_xvec,anly_struct.fAHP_Vm_mean,...
        struct('x','','y','(mV)','ttl',''),color_2,true,'$V_{\mathrm{fAHP,mean}}$');
    
    plot_with_smoothing(anly_xvec,anly_struct.AP_Vm_mean,...
        struct('x','','y','(mV)','ttl',''),color_1,true,'$V_{\mathrm{AP,mean}}$');
    set(gca,'xcolor', 'none','ycolor','k'); 
    
    legend_on_selected(gca);
    
    subplot(6,3,17); hold on;
   
    plot_with_smoothing(anly_xvec,anly_struct.Vthres_first,...
        struct('x','','y','(mV)','ttl',''),color_2,true,'$V_{\mathrm{thres}}^{(1)}$ (first thres)');
    
    plot_with_smoothing(anly_xvec,anly_struct.Vthres_mean_afterfirst,...
        struct('x','','y','(mV)','ttl',''),color_1,true,'$V_{\mathrm{thres,mean}}^{(> 1)}$ (mean after first)');
    set(gca,'xcolor', 'none','ycolor','k'); 
    
    legend_on_selected(gca);
        
    subplot(6,3,18); hold on;
   
    plot_with_smoothing(anly_xvec,anly_struct.Vbase,...
        struct('x','','y','(mV)','ttl',''),color_3,false,'$V_{\mathrm{pre-pulse}}$');
    
    plot_with_smoothing(anly_xvec,anly_struct.Vrest_1,...
        struct('x','','y','(mV)','ttl',''),color_2,false,'$V_{\mathrm{rest-1}}$ (Tuan)');
    
    plot_with_smoothing(anly_xvec,anly_struct.Vrest_2,...
        struct('x','','y','(mV)','ttl',''),color_1,false,'$V_{\mathrm{rest-2}}$ (Christian)');
    
    set(gca,'xcolor', 'none','ycolor','k'); 
    
    
    legend_on_selected(gca);
    
    pause(0.1); 
    export_fig(fullfile(fig_path, select_fignames{i}), '-dpng', '-r300', '-p0.02');
    close; 
    
    t1 = toc(t0)/60;
    fprintf('+ Cell "%s" done. Progress: %d/%d done. Elapsed %.1f min. ETA %.1f min.\n', ...
        select_fignames{i}, i, select_num_cells, t1, select_num_cells*t1/i - t1);
end

%% Customized functions for plotting

function plot_representative_traces(ax_list,repr_struct,stim_info,color_list,lgnd_names,repr_lbls)
time_offset = [-50, 75]; 
t_dep_rng = [stim_info.pulse_on_time,stim_info.pulse_off_time] + time_offset;
t_hyp_rng = [stim_info.hyperpol_on_time,stim_info.hyperpol_off_time] + time_offset;

t_ranges = {t_dep_rng, t_hyp_rng};

repr_fields = fieldnames(repr_struct);

for i = 1:length(t_ranges)
    
    hold(ax_list(i), 'on');
    set(ax_list(i), 'colororder', color_list);
    t_rng = t_ranges{i};
    
    for j = 1:length(repr_fields)
        field_repr = repr_fields{j};
        val_repr = repr_struct.(field_repr);
        t = val_repr.t;
        Vm = val_repr.Vm;
        
        filt_time = t >= t_rng(1) & t <= t_rng(2);
        plot(ax_list(i), t(filt_time)-t_rng(1), Vm(filt_time,:));
    end
    title(ax_list(i), repr_lbls.ttl{i}); 
    ylabel(ax_list(i), repr_lbls.y); 
    xlabel(ax_list(i), repr_lbls.x);
    xlim(ax_list(i), [0,diff(t_rng)]);
    
    if i == 1
        legend(ax_list(i), lgnd_names,'NumColumns',length(lgnd_names));
    end

end
end


function plot_with_smoothing(x,y,lbls,c,m_on,lgnd_name,sw,lsty)

if ~exist('c','var'), c=[0,0,0]; end
if ~exist('sw','var'), sw = 10; end
if ~exist('lsty','var'), lsty = '-'; end
if ~exist('m_on','var'), m_on = true; end
if ~exist('lgnd_name','var'), lgnd_name = ''; end
if ~exist('lbls','var'), lbls = struct('x','','y','','ttl',''); end

y_base = smoothdata(y(x<0), 'gaussian', sw);
y_post = smoothdata(y(x>0), 'gaussian', sw);


if ~isempty(lgnd_name)
    lgdn_tag = {'tag', 'legend-on', 'displayname', lgnd_name};
else
    lgdn_tag = {};
end

lineopts = {'color', c, 'linestyle', lsty};

plot(x(x<0),y_base,lineopts{:},lgdn_tag{:});
plot(x(x>0),y_post,lineopts{:});
xline(0, ':', 'linewidth', 2, 'color', 0.7*ones(1,3));

if m_on
    scatter(x,y,50,c(1:3),'filled', 'markerfacealpha',0.2,'markeredgecolor','none');
end

if ~isempty(lbls.x), xlabel(lbls.x); end
if ~isempty(lbls.y), ylabel(lbls.y); end
if ~isempty(lbls.ttl), title(lbls.ttl); end
set(gca, 'ycolor', c);

xlim(x([1,end]));
end

function legend_on_selected(ax)
if ~exist('ax','var'), ax = gca; end
legend(ax,findobj(ax,'tag', 'legend-on'),'FontSize',15);
end

