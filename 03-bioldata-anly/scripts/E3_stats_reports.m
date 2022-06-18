clc; clear; close all;
run start_up.m

%% Define paths
main_data_path = 'data/ip-data';
extra_data_path = 'data/extra-s1/processed/hansel-pham-extra-s1-vthres.mat';
v1_data_path = 'data/v1-data/v1-invivo_Li2020.mat';

fig_path = 'figures/ip-data/sum'; 
if ~exist(fig_path, 'dir')
    mkdir(fig_path); 
end

PROCESS_path = main_data_path; 
STATS_path = fullfile(main_data_path, 'stats');
if ~exist(STATS_path, 'dir')
    mkdir(STATS_path); 
end

ip_s1_data = load(fullfile(PROCESS_path, 'ip-select-data-analysis.mat'));
extra_s1_data = load(extra_data_path);

if check_if_has_v1(v1_data_path)
    li_v1_data = load(v1_data_path);
end

%% Define files
data_path = fullfile(STATS_path, 'ip-s1-data'); 
stat_path = fullfile(STATS_path, 'ip-s1-stat.xlsx'); 

num_sep = 4;

var_map = struct(...
    'ExpGroup', 'ExpGroup[S1]', ...
    'num_spikes_change', 'nspk_change[S1]', ...
    'dVthres_rest_change', 'dVTR_change[S1]', ...
    'Vthres_first_change', 'VT_change[S1]', ... 
    'Vrest_1_change',  'VR_change[S1]');

%% Write table of S1
s1_tbl = ip_s1_data.pooled_table;
writetable(s1_tbl, [data_path '.csv']);
writetable(s1_tbl, [data_path '.xlsx']);

%% Write table of V1
if check_if_has_v1(v1_data_path)
    v1_tbl = struct2table(structfun(@(x) x.data, li_v1_data.data, 'uni', 0));
    writetable(v1_tbl, fullfile(STATS_path, 'v1-data.csv'));
    writetable(v1_tbl, fullfile(STATS_path, 'v1-data.xlsx'));
else
    pseudo_v1 = struct(...
        'Vrest', nan(10,1), ...
        'Vthres', nan(10,1), ...
        'dV', nan(10,1));
    v1_tbl = struct2table(pseudo_v1); 
end

%% Correlation test between response variables
g_name = 'ExpGroup';
Y_names = fieldnames(var_map);
Y_names = Y_names(~contains(Y_names, g_name));

response_data = cell2mat(cellfun(@(x) s1_tbl.(x), Y_names', 'UniformOutput', false));
response_vars = cellfun(@(x) var_map.(x), Y_names, 'UniformOutput', false);
group_data = s1_tbl.(g_name); 
unq_groups = [{'all'}; unique(group_data)];

m2T = @(m, name, group) array2table(m, ...
    'RowNames', response_vars, ...
    'VariableNames', response_vars, ...
    'DimensionNames', {sprintf('%s [group=%s]', name, upper(group)), 'Variables'});

output_common_args = {stat_path, 'WriteRowNames', true, 'Sheet', 'corr-nspk-and-Vm'};

curr_row = 1; 

for i = 1:length(unq_groups)
    grp_i= unq_groups{i};
    if strcmp(grp_i, 'all')
        loc_i = true(size(response_data,1), 1);
    else
        loc_i = strcmp(group_data, grp_i);
    end
    
    [r, p] = corrcoef(response_data(loc_i, :));
    T_rho = m2T(r, 'rho', grp_i);
    T_pval = m2T(p, 'pval', grp_i);
    
    curr_col = 'A';
    writetable(T_rho, output_common_args{:}, 'Range', sprintf('%c%d', curr_col, curr_row));
    
    curr_col = curr_col + width(T_rho) + num_sep;
    writetable(T_pval, output_common_args{:}, 'Range', sprintf('%c%d', curr_col, curr_row));
    
    curr_row = curr_row + height(T_rho) + num_sep + 2; 

end

%% T-Tests pre-vs-post within each groups separately
response_roots = cellfun(@(y) strrep(y, '_change', ''), Y_names, 'uni', 0);
num_res_vars = length(response_roots);

output_common_args = {stat_path, 'WriteRowNames', true, 'Sheet', 'ttest-pre-vs-post'};

curr_row = 1;
for i = 1:length(unq_groups)
    grp_i = unq_groups{i};
    if strcmp(grp_i, 'all')
        loc_i = true(size(response_data,1), 1);
    else
        loc_i = strcmp(group_data, grp_i);
    end
    stat_res = cell(num_res_vars, 1);
    for j = 1:length(response_roots)
        field_j = response_roots{j};
        base_data = get_base(s1_tbl, field_j);
        post_data = get_post(s1_tbl, field_j);
        pair_data = {base_data(loc_i), post_data(loc_i)};
        stat_res{j} = perform_stat_tests(pair_data{:});
        diff_data = pair_data{2} - pair_data{1};
         
        stat_res{j}.n = length(diff_data);
        stat_res{j}.mean_diff = mean(diff_data);
        stat_res{j}.std_diff = std(diff_data);
    end
    stat_res = struct2table(vertcat(stat_res{:}), ...
        'RowNames', cellfun(@(x) strrep(var_map.([x '_change']), '_change', ''), response_roots, 'uni', 0), ... lazy, dont judge
        'DimensionNames', {sprintf('T-test pre-vs-post [group=%s]', upper(grp_i)), 'Variables'});
    
    writetable(stat_res, output_common_args{:}, 'Range', sprintf('A%d', curr_row));
    curr_row = curr_row + height(stat_res) + num_sep; 
end

%% ANOVA-1 
output_common_args = {stat_path, 'Sheet', 'anova1-Grp-vs-Vm-nspk'};
use_sign_opts = [false, true]; 
cut_off_sign = 0.01; 

for k = 1:length(use_sign_opts)
    use_sign = use_sign_opts(k);
    
    sheet_name = output_common_args{3};
    if use_sign
        sheet_name = sprintf('%s_use-sign', sheet_name);
    else
        sheet_name = sprintf('%s_use-vals', sheet_name);
    end
    output_common_args_k = {output_common_args{1:2}, sheet_name};
    
    curr_row = 1;
    
    for i = 1:length(Y_names)
        y_name = Y_names{i};
        y_vals = s1_tbl.(y_name); 
        
        if use_sign
            y_vals = sign(y_vals) .* (abs(y_vals) > cut_off_sign);
        end
        
        [~, T_anova1, stats] = anova1(y_vals, s1_tbl.(g_name), 'off');
        T_anova1{1,1} = sprintf('ANOVA-1, Y = %s', var_map.(y_name));
        T_anova1{2,1} = var_map.(g_name);
        T_anova1 = cell2table_withnames(T_anova1);
        
        [c,~,~,gnames] = multcompare(stats, 'display', 'off');
        pair_ids = arrayfun(@(x) sprintf('pair-%d', x), 1:size(c,1), 'uni', 0)';
        mcmp_tbl = cell2table([pair_ids, gnames(c(:,1)), gnames(c(:,2)), num2cell(c(:,3:6))], ...
            'VariableNames', {'Compare', 'Group-1', 'Group-2', 'low[ci95]_diff', 'mean_diff', 'upp[ci95]_diff', 'p[diff=0]'});
        
        writetable(T_anova1, output_common_args_k{:}, 'Range', sprintf('A%d', curr_row));
        writetable(mcmp_tbl, output_common_args_k{:}, 'Range', sprintf('H%d', curr_row));
        
        curr_row = curr_row + max([height(T_anova1), height(mcmp_tbl)]) + num_sep;
    end
end

%% MANOVA-1

[d_manova,p_manova] = manova1(response_data, group_data);
row_names = [{'input: X', 'input: Group', 'output: d'}, arrayfun(@(i) sprintf('output: p[d=%d]', i-1), 1:length(p_manova), 'uni', 0)]'; 
x_desc = ['S1 changes in ' strjoin( regexprep(response_vars,'(_change|\[S1\])',''), ', ')];
table_vals = [{x_desc; 'ExpGroup[S1]'}; {d_manova}; mat2colcell(p_manova)];
T_manova = table(row_names, table_vals, 'VariableNames', {'MANOVA-1', 'Value'});

writetable(T_manova, stat_path, 'Sheet', 'manova1-Grp-vs-Vm-nspk');

%% ANOVA-2

output_common_args = {stat_path, 'Sheet', 'anova2-nspk-vs-Vm-and-Grp'};
[~, T_anovan] = anovan(s1_tbl.num_spikes_change, ...
    {s1_tbl.ExpGroup, s1_tbl.Vrest_1_change, s1_tbl.Vthres_first_change}, ...
    'varnames', {'ExpGroup[S1]','VR_change[S1]', 'VT_change[S1]'}, ...
    'continuous', [2,3], 'model',2, 'display', 'off');
T_anovan{1,1} = sprintf('ANOVA-2, Y = %s', var_map.num_spikes_change);
writetable(cell2table_withnames(T_anovan), output_common_args{:});

curr_row = size(T_anovan,1);

[~, T_anovan] = anovan(s1_tbl.num_spikes_change, ...
    {s1_tbl.ExpGroup, s1_tbl.dVthres_rest_change}, ...
    'varnames', {'ExpGroup[S1]','dVTR_change[S1]'}, ...
    'continuous', 2, 'model', 2, 'display', 'off');
T_anovan{1,1} = sprintf('ANOVA-2, Y = %s', var_map.num_spikes_change);
writetable(cell2table_withnames(T_anovan), output_common_args{:}, ...
    'Range', sprintf('A%d', curr_row + num_sep));

%% Write summary stats and comparing pre/post or between S1/V1 to a latex table
if ~check_if_has_v1(v1_data_path)
    warning('For V1 summary stats and tests with S1, will use a pseudo data table in order to finish the table construction');
end

data_source_mapping = struct;
data_source_mapping.V_R = struct('V1', 'Vrest', 'S1', 'Vrest_1');
data_source_mapping.V_T = struct('V1', 'Vthres', 'S1', 'Vthres_first');
data_source_mapping.dV_TR = struct('V1', 'dV', 'S1', 'dVthres_rest');

data_source_latex = struct;
data_source_latex.V_R = '$V_{\mathrm{R}}$';
data_source_latex.V_T = '$V_{\mathrm{T}}$';
data_source_latex.dV_TR = '$\Delta V_{\mathrm{TR}}$';

latex_table_header = {...
    '\documentclass{article}', ...
    '\usepackage{multirow,caption,float,fontsize}', ...
    '\usepackage[a4paper,margin=1in,landscape]{geometry}',...
    '\renewcommand{\arraystretch}{1.25}', ...
    '\captionsetup{aboveskip=5pt,labelfont=bf,justification=justified,labelsep=period}', ...
    '\changefontsize[15pt]{14pt}', ...
    '\begin{document}', ...
    '\begin{table}[H]', ...
    '\centering', ...
    '\caption{S1 \textit{in vitro} (with \textbf{IP}) and V1 \textit{in vivo} data}', ...
    '\begin{tabular}{|l|l|l|l|l|l|l|l|l|l|}',...
    '\hline Variable & \multicolumn{2}{l}{Source} & range & mean $\pm$ SD & $n$ & $t$-stat & $p_t$ & $F$-stat & $p_F$ \\', ...
    '\hline'
};

latex_table_footer = {...
    '\end{tabular}', ...
    '\end{table}',  ...
    '\end{document}' 
};

common_args = {s1_tbl, v1_tbl, data_source_mapping, data_source_latex};
data_latex_rows = cellfun(@(x) generate_latex_rows(x, common_args{:}), fieldnames(data_source_mapping), 'uni', 0);
data_latex_rows = horzcat(data_latex_rows{:});
latex_table = [latex_table_header, data_latex_rows, latex_table_footer];

fileID = fopen(fullfile(STATS_path, 'stats-table-s1-v1.tex'),'w');
for i = 1:length(latex_table)
    row = latex_table{i};
    fprintf(fileID, '%s\n', row);
end
fclose(fileID);

%% Helper functions

function T = cell2table_withnames(C, include_rownames)
    row_args = {};
    if ~exist('include_rownames', 'var') 
        include_rownames = false; 
    end
    if include_rownames 
        row_args = {'RowNames', C(2:end,1)};
    end
    
    start_ind = 1 + include_rownames;
    
    T = cell2table(C(2:end,start_ind:end), ...
        'VariableNames', C(1,start_ind:end), ...
        row_args{:}); 
end

function latex_rows = generate_latex_rows(data_source, s1_tbl, v1_tbl, data_source_mapping, data_source_latex)
    data_latex_name = data_source_latex.(data_source);

    pairs_to_test = struct(...
        's1_pre_vs_post', {{'s1_base', 's1_post'}}, ...
        's1_vs_v1', {{'s1_all', 'v1'}}...
    );

    v1_field = data_source_mapping.(data_source).V1;
    s1_field = data_source_mapping.(data_source).S1;

    data = struct(...
        's1_base', get_base(s1_tbl, s1_field), ...
        's1_post', get_post(s1_tbl, s1_field), ...
        's1_all', get_base_and_post(s1_tbl, s1_field), ...
        'v1', v1_tbl.(v1_field) ...
    );

    data_summary = structfun(@get_summary_stats, data, 'uni', 0);
    data_stattests = structfun(@(x) perform_stat_tests(data.(x{1}), data.(x{2})), pairs_to_test, 'uni', 0);

    latex_row_1 = [...
        {tex_multirow(4, data_latex_name), tex_multirow(3, 'S1'), 'pre'}, ...
        tex_summary(data_summary.s1_base), ...
        tex_stattest(data_stattests.s1_pre_vs_post), ...
        ];

    latex_row_2 = [...
        {' ', ' ', 'post'}, ...
        tex_summary(data_summary.s1_post), ...
        tex_stattest([]), ...
        ];

    latex_row_3 = [...
        {' ', ' ', 'all'}, ...
        tex_summary(data_summary.s1_all), ...
        tex_stattest(data_stattests.s1_vs_v1), ...
        ];

    latex_row_4 = [...
        {' ', tex_multicol(2, 'V1')}, ...
        tex_summary(data_summary.v1), ...
        tex_stattest([]), ...
        ];

    latex_rows = {...
        tex_finalize_row(latex_row_1), ...
        tex_finalize_row(latex_row_2), ...
        tex_finalize_row(latex_row_3, '\cline{7-10}'), ...
        tex_finalize_row(latex_row_4, '\cline{2-3}', '\hline')};
end

function r = tex_finalize_row(r, prefix, suffix)
    if ~exist('prefix', 'var'), prefix = ''; end
    if ~exist('suffix', 'var'), suffix = ''; end
    r = join(r, '&');
    r = [prefix, r{1}, '\\', suffix];
end

function r = tex_multirow(n, s, a)
    if ~exist('a', 'var'), a = '*'; end
    r = sprintf('\\multirow{%d}{%s}{%s}', n, a, s); 
end

function r = tex_multicol(n, s, a)
    if ~exist('a', 'var'), a = 'l|'; end
    r = sprintf('\\multicolumn{%d}{%s}{%s}', n, a, s); 
end

function r = tex_bounds(data)
    bounds = data.bounds;
    r = sprintf('$[%.1f, %.1f]$', bounds(1), bounds(2));
end

function r = tex_mean(data)
    r = sprintf('$%.2f \\pm  %.2f$', data.mean, data.std);
end

function r = tex_sample(data)
    r = sprintf('$%d$', data.n);
end

function r = tex_summary(data)
    r = {...
        tex_bounds(data), ...
        tex_mean(data), ...
        tex_sample(data)
    };
end

function r = tex_pval(p, alpha)
    if ~exist('alpha', 'var'), alpha = 0.05; end
    if p < alpha 
        if p >= 1e-4
            p_txt = sprintf('\\mathbf{%.4f}', p);
        else 
            p_txt = sprintf('\\mathbf{%.5g}', p);
        end
    else
        p_txt = sprintf('%.4f', p);
    end

    star_levels = [0.05, 0.01, 1e-3]; 
    signif_level = repmat('*', [sum(p < star_levels), 1]);

    if ~isempty(signif_level)
        r = sprintf('$%s^{%s}$', p_txt, signif_level);
    else
        r = sprintf('$%s$', p_txt);
    end
end

function r = tex_statval(s)
    r = sprintf('$%.4f$', s);
end 

function r = tex_stattest(data)
    if isempty(data)
        r = {' ', ' ', ' ', ' '}; 
        return 
    end

    r = {...
        tex_multirow(2, tex_statval(data.t_stat)), ...
        tex_multirow(2, tex_pval(data.p_t)), ...
        tex_multirow(2, tex_statval(data.f_stat)) ...
        tex_multirow(2, tex_pval(data.p_f)), ...
    };

end

function r = get_summary_stats(x)
    r = struct;
    r.n = length(x); 
    r.bounds = [min(x), max(x)];
    r.median = median(x); 
    r.mean = mean(x); 
    r.std = std(x); 
    r.sem = r.std / sqrt(r.n);
end

function r = perform_stat_tests(x1, x2)
    r = struct;

    ttest_fn = @ttest;
    if length(x1) ~= length(x2)
        ttest_fn = @ttest2;
    end
    [~, r.p_t, ~, tstats] = ttest_fn(x1, x2); 
    r.t_stat = tstats.tstat;

    [~, r.p_f, ~, fstats] = vartest2(x1, x2);
    r.f_stat = fstats.fstat; 
end

function v = get_base(T,f)
    v = T.([f '_base']);
end

function v = get_post(T,f)
    v = T.([f '_post']);
end

function v = get_base_and_post(T,f)
    v = [get_base(T,f);get_post(T,f)];
end

function has_v1_data = check_if_has_v1(data_path)
    has_v1_data = exist(data_path, 'file') ~= 0;
    if ~has_v1_data
        warning('Does not have "%s" for V1 data, will skill this', data_path);
    end
end