function res = process_stats(X, d, sq)
% take mean and sem assuming the last dimension is # samples
% unless otherwise specified by using `d`
% d = -1 is equivalent of last dim
% sq is to decide whether to further squeeze, def = false 

    sz = size(X);
    if ~exist('sq', 'var'), sq = false; end 
    if ~exist('d', 'var'), d = -1; end
    if d==-1, d = length(sz); end 
    
    n = size(X,d); 
    sqrt_n = sqrt(n);
    res = struct(...
        'mean', mean(X,d),...
        'sem', std(X,[],d) / sqrt_n);
    
    if sq
        res = structfun(@squeeze, res, 'uni', 0);
    end
end