clc; clear; close all;
run start_up.m

extra_s1_data = load('data/extra-s1/processed/hansel-pham-extra-s1-vthres.mat');
gill_s1_data = load('data/gill-data/Layer II.III S1 BC IP/processed/gill-select-data-analysis.mat');

%% get summary of s1
s1_tbl = gill_s1_data.pooled_table;
summary(s1_tbl) ;

%% t test for pre-vs-post


[h,p,ci,stats] = ttest(s1_tbl.Vrest_1_base, s1_tbl.Vrest_1_post);
[h,p,ci,stats] = ttest(s1_tbl.Vthres_first_base, s1_tbl.Vthres_first_post);
[h,p,ci,stats] = ttest(s1_tbl.dVthres_rest_base, s1_tbl.dVthres_rest_post);

%% f test for pre-vs-post 
[h,p,ci,stats] = vartest2(s1_tbl.Vrest_1_base, s1_tbl.Vrest_1_post);
[h,p,ci,stats] = vartest2(s1_tbl.Vthres_first_base, s1_tbl.Vthres_first_post);
[h,p,ci,stats] = vartest2(s1_tbl.dVthres_rest_base, s1_tbl.dVthres_rest_post);

%% v1 data 
li_v1_data = load('~/Documents/MATLAB/MATLABDrive/HanselLab/matlab/code/combo_and_thres_off/data/V1invivo_Li2020.mat');
v1_tbl = struct2table(structfun(@(x) x.data, li_v1_data.data, 'uni', 0));

[~,p,~,stats]=ttest2([s1_tbl.Vrest_1_base;s1_tbl.Vrest_1_post], v1_tbl.Vrest);
[~,p,~,stats]=ttest2([s1_tbl.Vthres_first_base;s1_tbl.Vthres_first_post], v1_tbl.Vthres);
[~,p,~,stats]=ttest2([s1_tbl.dVthres_rest_base;s1_tbl.dVthres_rest_post], v1_tbl.dV);

%%
unq_groups = unique(s1_tbl.ExpGroup);

x_var_names = {'Vrest_1_change', 'Vthres_first_change', 'dVthres_rest_change'};
y_var_name = 'num_spikes_change';

for i = 1:length(x_var_names)
    x_var_name = x_var_names{i};
    x_vec = s1_tbl.(x_var_name);
    y_vec = s1_tbl.(y_var_name);
    
    fprintf('CORR TEST: "%s" vs "%s"\n', x_var_name, y_var_name);
    
    for j = 1:length(unq_groups)
        grp_j = unq_groups{j};
        loc_j = strcmp(s1_tbl.ExpGroup, grp_j);
        
        [r,p] = corrcoef(x_vec(loc_j), y_vec(loc_j));
        r_val = r(1,2);
        p_val = p(1,2); 
        
        if p_val < 0.05 
            p_fmt = '<strong>%.3g</strong>';
        else
            p_fmt = '%.3g';
        end
        
        fprintf(['\t + GROUP ["%s"]: (unadj) p-value = ' p_fmt ', rho = %.2f \n'], grp_j, p_val, r_val);
    end
end

%%

s1_tbl.exp_elec = repmat("none", [size(s1_tbl,1), 1]);
s1_tbl.exp_elec(contains(s1_tbl.ExpGroup,'elec')) = "elec";
s1_tbl.exp_elec = categorical(s1_tbl.exp_elec);

s1_tbl.exp_chol = repmat("none", [size(s1_tbl,1), 1]);
s1_tbl.exp_chol(contains(s1_tbl.ExpGroup,'chol')) = "chol";
s1_tbl.exp_chol = categorical(s1_tbl.exp_chol);

%%
anovan(s1_tbl.num_spikes_change, ...
    {s1_tbl.exp_chol, s1_tbl.exp_elec, ...
    s1_tbl.Vrest_1_change, s1_tbl.Vthres_first_change}, ...
    'continuous', [3,4], ...
    'varnames',{'IsChol','IsElec','V_R change', 'V_T change'}, ...
    'model',1)

% anovan(s1_tbl.num_spikes_change, ...
%     {s1_tbl.exp_chol, s1_tbl.exp_elec, ...
%     s1_tbl.Vrest_1_change, s1_tbl.Vthres_first_change}, ...
%     'continuous', [1,2,3,4], ...
%     'varnames',{'IsChol','IsElec','V_R change', 'V_T change'}, ...
%     'model',2)

%% anovan need a control group to do this correctly
% s1_tbl.exp_elec = contains(s1_tbl.ExpGroup,'elec');
% s1_tbl.exp_chol = contains(s1_tbl.ExpGroup,'chol');
% 
% s1_tbl.exp_elec = s1_tbl.exp_elec + 1e-10 * rand(size(s1_tbl.exp_elec));
% s1_tbl.exp_chol = s1_tbl.exp_chol + 1e-10 * rand(size(s1_tbl.exp_chol));
% 
% anovan(s1_tbl.num_spikes_change, ...
%     {s1_tbl.exp_chol, s1_tbl.exp_elec, ...
%     s1_tbl.Vrest_1_change, s1_tbl.Vthres_first_change}, ...
%     'continuous', [1,2,3,4], ...
%     'varnames',{'IsChol','IsElec','V_R change', 'V_T change'}, ...
%     'model',1)

%%

anovan(s1_tbl.num_spikes_change, ...
    {s1_tbl.ExpGroup, ...
    s1_tbl.Vrest_1_change, s1_tbl.Vthres_first_change}, ...
    'continuous', [2,3], ...
    'varnames',{'ExpGroup','V_R change', 'V_T change'}, ...
    'model',1)

%%

s1_tbl.exp_chol = contains(s1_tbl.ExpGroup,'chol');

anovan(s1_tbl.num_spikes_change, ...
    {s1_tbl.exp_chol, ...
    s1_tbl.Vrest_1_change, s1_tbl.Vthres_first_change}, ...
    'continuous', [2,3], ...
    'varnames',{'IsChol','V_R change', 'V_T change'}, ...
    'model',2)

%%
s1_tbl.exp_elec = contains(s1_tbl.ExpGroup,'elec');

anovan(s1_tbl.num_spikes_change, ...
    {s1_tbl.exp_chol, ...
    s1_tbl.Vrest_1_change, s1_tbl.Vthres_first_change}, ...
    'continuous', [2,3], ...
    'varnames',{'IsElec','V_R change', 'V_T change'}, ...
    'model',2)

