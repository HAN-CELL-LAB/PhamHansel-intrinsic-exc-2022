function [h,p,chi2,df] = chi2tod(X, alpha)
    % chi2 test of independence 
    % X = [90,60,104,95; 30,50,51,20;30,40,45,35]
    if ~exist('alpha', 'var')
        alpha = 0.05; 
    end
    
    [nrow, ncol] = size(X);
    df = (nrow - 1) * (ncol - 1); 
    
    E = sum(X,1) .* sum(X,2) / sum(X, 'all'); 
    chi2 = sum((X - E).^2 ./ E, 'all');
    p = 1 - chi2cdf(chi2,df);
    h = p < alpha;
end
