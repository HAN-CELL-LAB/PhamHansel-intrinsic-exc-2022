function c = count_from_unique_id(x, max_x)
c = zeros(max_x+1, 1); 
[unq_x,~,icx] = unique(x);
c(unq_x+1) = accumarray(icx,1);
end