function [p, k] = id2prob(v) 
    p = groupcounts(v); 
    p = p / sum(p); 
    k = length(p); 
end
