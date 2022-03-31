function [Ybid_unqcnt, pattern_ids, I] = combo_and_result(N_I,N_X,genWfun,num_runs,num_thres)

pattern_ids = 0:(2^N_I-1);
n_p = length(pattern_ids); 
I = de2bi(pattern_ids);

thres_X = linspace(0,1,num_thres);
thres_Y = linspace(0,3,num_thres);

Ybid_unqcnt = cell(num_runs, 1);

% scale_W_YX = N_I/N_X;
scale_W_YX = 1;
for n = 1:num_runs
    W_XI = genWfun(N_X,N_I);
    W_XI = W_XI ./ sum(W_XI, 2);
    W_YX = scale_W_YX * W_XI';
    
    Xpre = W_XI * I';
    X = Xpre > reshape(thres_X, [1,1,num_thres]);
    
    Y = cellfun(@(x) W_YX * x, mat2cell_splitdim(X,3), 'uni', 0);
    Y = cat(3, Y{:}) > reshape(thres_Y, [1,1,1,num_thres]);
    
    Y = squeeze(mat2cell(Y, size(Y,1), size(Y,2), ones(size(Y,3),1), ones(size(Y,4),1))); 
 
    Yb = cellfun(@(y) accumarray(bi2de(y')+1,1,[n_p,1])/n_p, Y, 'uni', 0);
  
    Ybid_unqcnt{n} = cell2mat(reshape(Yb, [1, size(Yb)]));
    
end

dim_Ybidunqcnt = length(size(Ybid_unqcnt{1})); 
cat_dim = dim_Ybidunqcnt+1;
Ybid_unqcnt = mean(cat(cat_dim,Ybid_unqcnt{:}), cat_dim);

end