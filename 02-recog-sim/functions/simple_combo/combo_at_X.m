function nYact_per_k = combo_at_X(N,genWfun,num_runs,num_thres) 
W = genWfun(1,N); 
W = W/sum(W);

thres = linspace(0,1,num_thres); 

X = zeros(N,num_runs);
for i = 1:num_runs
    k_X = randperm(N+1,1)-1; 
    X(randperm(N,k_X),i) = 1; 
end

Y = W * X;
Y = Y' > thres;

k_X = sum(X, 1);

unq_k = 0:N;
nYact_per_k = arrayfun(@(k) mean(Y(k_X == k,:),1), unq_k, 'uni', 0);
nYact_per_k = vertcat(nYact_per_k{:});

% res = struct; 
% res.unq_k = unq_k; 
% res.nYact = nYact; 
% res.nYact_per_k = nYact_per_k; 
end