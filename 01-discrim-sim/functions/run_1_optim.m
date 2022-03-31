function results = run_1_optim(N_X, N_Y, n_iter, eta_grad, epsi_deltafun)
    if ~exist('n_iter', 'var'), n_iter = 2500; end
    if ~exist('eta_grad', 'var'), eta_grad = 1e0; end
    if ~exist('epsi_deltafun', 'var'), epsi_deltafun = 1e-2; end


    % initialize results struct
    [X_data, ~, Y_data, ~] = assign_XY_pairs(N_X, N_Y);

    results = struct(...
        'optim_progress', struct(...
            'loss', zeros(n_iter, 1), ...
            'entropy', zeros(n_iter, 1), ...
            'num_unq', zeros(n_iter, 1)), ...    
        'net_config', struct(...
            'init', [], ...
            'best', []) ...
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

    results = results.net_config;
end

function net_config = save_netconfig(net_config, iter, loss, ent, unq)
    net_config.iter = iter;
    net_config.loss = loss;
    net_config.entropy = ent; 
    net_config.num_unq = unq;
end
