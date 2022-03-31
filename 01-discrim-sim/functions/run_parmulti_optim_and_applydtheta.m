function results = run_parmulti_optim_and_applydtheta(saved_filename, n_sim, N_X, N_Y, dtheta_vec, varargin)
    % TODO: add num_tryouts to args
    num_tryouts = 5;
    
    configurations = struct(...
        'n_sim', n_sim, ...
        'N_X', N_X, ...
        'N_Y', N_Y, ...
        'dtheta_vec', dtheta_vec);
    
    simulations = cell(n_sim, 1);
    parfor i_sim = 1:n_sim
        
        i_try = 0;
        success_ith = false; 
        
        while ~success_ith && i_try < num_tryouts
            try
                simulations{i_sim} = run_1_optim_and_applydtheta(N_X, N_Y, dtheta_vec, varargin{:});
                success_ith = true;
            catch
                simulations{i_sim} = [];
            end
            i_try = i_try + 1;            
        end
    end
    sucess_sims = cellfun(@(x) ~isempty(x), simulations);
    simulations = simulations(sucess_sims);
    
    fprintf('(%03d / %03d sims) \t ',  sum(sucess_sims), n_sim);
    
    L1_fields = fieldnames(simulations{1}.apply_dtheta);
    L2_fields = fieldnames(simulations{1}.apply_dtheta.(L1_fields{1}));
    template_L2 = empty_struct_like(L2_fields); 
    
    analyses = empty_struct_like(L1_fields, template_L2); 
    
    for i1 = 1:length(L1_fields)
        F1 = L1_fields{i1};
        for i2 = 1:length(L2_fields)
            F2 = L2_fields{i2};
            
            mat_aggr = cellfun(@(x) x.apply_dtheta.(F1).(F2), simulations, 'uni', 0);
            tmp_val = mat_aggr{1}(1);
            
            if isnumeric(tmp_val) 
                mat_aggr = horzcat(mat_aggr{:});
            elseif iscell(tmp_val)
                if ~strcmpi(F2, 'prob')
                    warning('Currently only process "%s" for cell type. Will ignore', F2);
                    continue; 
                end
                
                max_len = max(cellfun(@(x) max(cellfun(@length, x)), mat_aggr, 'uni', 1));
                mat_aggr = cellfun(@(x) process_probcells(x, max_len), mat_aggr, 'uni', 0);
                mat_aggr = cat(3, mat_aggr{:});
            else
                warning('Type of values in field "%s->%s" cannot proceed. Will ignore', F1, F2);
                continue;
            end
            
            analyses.(F1).(F2) = process_stats(mat_aggr);
        end
    end
    
    save(saved_filename, 'configurations', 'simulations', 'analyses');
    results = analyses;         
end

function m = process_probcells(m, L)
m = cellfun(@(x) to_row_vec(extendvec_fillmissing(x, L)), m, 'uni', 0);
m = vertcat(m{:});
end