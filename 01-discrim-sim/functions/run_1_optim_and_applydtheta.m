function results = run_1_optim_and_applydtheta(N_X, N_Y, dtheta_vec, n_iter, eta_grad, epsi_deltafun)
    if ~exist('n_iter', 'var'), n_iter = 5000; end
    if ~exist('eta_grad', 'var'), eta_grad = 1e0; end
    if ~exist('epsi_deltafun', 'var'), epsi_deltafun = 1e-2; end

    n_dtheta = length(dtheta_vec); 

    % initialize results struct
    [X_data, X_lbls, Y_data, Y_lbls] = assign_XY_pairs(N_X, N_Y);

    template_applydtheta_results = struct(...
        'prob', {cell(n_dtheta, 1)}, ...
        'entropy', zeros(n_dtheta, 1), ...
        'num_unq', zeros(n_dtheta, 1));

    results = struct(...
        'sim_config', struct(...
            'N_X', N_X, ...
            'N_Y', N_Y, ...
            'n_iter', n_iter, ...
            'eta_grad', eta_grad, ...
            'epsi_deltafun', epsi_deltafun, ...        
            'dtheta_vec', dtheta_vec, ...
            'X_labels', X_lbls, ...
            'Y_labels', Y_lbls), ...
        'optim_progress', struct(...
            'loss', zeros(n_iter, 1), ...
            'entropy', zeros(n_iter, 1), ...
            'num_unq', zeros(n_iter, 1)), ...    
        'net_config', struct(...
            'init', [], ...
            'best', []), ...
        'apply_dtheta', struct(...
            'init', template_applydtheta_results, ...
            'best', template_applydtheta_results) ...
        );  

    % init for optim
    W = rand(N_Y,N_X);
    W = W./sum(W,2);
    theta = rand(N_Y,1); 

    results.net_config.init = struct('W', W, 'theta', theta);
    results.net_config.best = struct('W', W, 'theta', theta);

    % optim
    for i_iter = 1:n_iter

        % infer
        pre_Y = W * X_data - theta;
        Y_infer = heavisidefun(pre_Y);

        Y_inf_ids = bi2de(Y_infer');
        err_Y = (Y_infer - Y_data); 

        % analyze
        curr_loss = mean(err_Y.^2, 'all');
        [p,k] = id2prob(Y_inf_ids);
        curr_ent = entropy(p);
        curr_unq = k;

        % save to vec
        results.optim_progress.loss(i_iter) = curr_loss;
        results.optim_progress.entropy(i_iter) = curr_ent;
        results.optim_progress.num_unq(i_iter) = curr_unq;

        % save to struct
        if i_iter == 1 
            results.net_config.init = save_netconfig(results.net_config.init, ...
                i_iter, curr_loss, curr_ent, curr_unq);
            results.net_config.best = save_netconfig(results.net_config.best, ...
                i_iter, curr_loss, curr_ent, curr_unq);
        else
            if results.net_config.best.entropy < curr_ent
                results.net_config.best = save_netconfig(results.net_config.best, ...
                    i_iter, curr_loss, curr_ent, curr_unq);

                results.net_config.best.W = W;
                results.net_config.best.theta = theta;
            end
        end

        % update 
        delta_preY = deltafun(pre_Y, epsi_deltafun);
        dL_dW = 2*(err_Y.*delta_preY)*X_data' / size(X_data,2);
        dL_dtheta = -2*mean(err_Y.*delta_preY, 2);

        W = W - eta_grad * dL_dW;
        theta = theta - eta_grad * dL_dtheta;
    end

    % apply dtheta

    net_fields = fieldnames(results.net_config);

    for i_field = 1:length(net_fields)
        field_name = net_fields{i_field};

        W = results.net_config.(field_name).W;
        theta = results.net_config.(field_name).theta;

        for i_dtheta = 1:n_dtheta
            % dtheta -> infer
            dtheta = dtheta_vec(i_dtheta);
            theta_new = theta + dtheta;
            Y_infer_new = heavisidefun(W * X_data - theta_new);
            Y_inf_ids_new = bi2de(Y_infer_new');

            % analyze
            [p,k] = id2prob(Y_inf_ids_new);
            results.apply_dtheta.(field_name).prob{i_dtheta} = sort(p, 'descend');
            results.apply_dtheta.(field_name).entropy(i_dtheta) = entropy(p);
            results.apply_dtheta.(field_name).num_unq(i_dtheta) = k;
        end

    end

end

function net_config = save_netconfig(net_config, iter, loss, ent, unq)
    net_config.iter = iter;
    net_config.loss = loss;
    net_config.entropy = ent; 
    net_config.num_unq = unq;
end
