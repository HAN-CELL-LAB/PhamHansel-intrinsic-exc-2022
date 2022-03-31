clc; clear; close all;
run start_up;

full_file_name = 'data/V1invivo_Li2020.csv';
tbl = readtable(full_file_name, 'CommentStyle', '/*');
tbl.dV_mV = tbl.Vthres_mV - tbl.Vrest_mV;
data = table2struct(tbl);
data = structarray_to_struct(data);

latex_convert = @(s) sprintf('%s_{\\mathrm{%s}}', s(1), s(2:end));
latex_unit = @(s) regexprep(s, {'Ohm'}, {'\\Omega'});
for fn = fieldnames(data)'
    field_name = fn{:};
    var_name_and_unit = strsplit(field_name, '_');
    var_name = var_name_and_unit{1};
    var_unit = latex_unit(var_name_and_unit{2});
    var_latex = latex_convert(var_name);
    var_data = data.(field_name);
    
    data = rmfield(data, field_name);
    data.(var_name) = struct(...
        'unit', var_unit, ...
        'latex', var_latex, ...
        'data', var_data);
    
end

data.dV.latex = '\Delta\mathrm{V}';
tbl_vars = tbl.Properties.VariableNames;
tbl_vals = tbl.Variables;

summary_stats = table(...
    mean(tbl_vals,1)', ...
    median(tbl_vals,1)', ...
    std(tbl_vals,1)', ...
    std(tbl_vals,1)'/sqrt(size(tbl_vals,1)), ...
    min(tbl_vals,[],1)', ...
    max(tbl_vals,[],1)', ...`
    'RowNames', tbl_vars, ...
    'VariableNames', {'mean', 'median', 'sd', 'sem', 'min', 'max'});

[file_path, file_name, ~] = fileparts(full_file_name);
save(fullfile(file_path, [file_name '.mat']), 'data', 'summary_stats');


array2table(table2array(summary_stats)', 'VariableNames', summary_stats.Row, 'RowNames', summary_stats.Properties.VariableNames)
