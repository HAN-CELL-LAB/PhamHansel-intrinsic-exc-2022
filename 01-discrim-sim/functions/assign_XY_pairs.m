function [X_data, X_lbls, Y_data, Y_lbls] = assign_XY_pairs(N_X, N_Y)
    num_Xs = 2^N_X;
    max_Ys = 2^N_Y - 1;
    X_lbls = (0:(num_Xs-1))';

    if N_Y < N_X
        Y_lbls = randi(max_Ys, num_Xs, 1);
    elseif N_X == N_Y 
        Y_lbls = randperm(num_Xs)' - 1;
    else 
        Y_lbls = randperm(max_Ys+1, num_Xs) - 1;
    end

    X_data = de2bi(X_lbls)';
    Y_data = de2bi(Y_lbls)';

end
