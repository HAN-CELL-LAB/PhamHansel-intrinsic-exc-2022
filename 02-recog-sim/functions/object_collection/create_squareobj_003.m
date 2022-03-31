function X = create_squareobj_003(N,k)
X = zeros(N);
if mod(N,2) ~= 0
    error('N must be even'); 
end

mid_X = numel(X)/2;
X(1:mid_X) = 1; 

if k == 0 
   return; 
end 

k = bound_minmax(k,1,mid_X);
swap_X = (mid_X+1-k):(mid_X+k); 
X(swap_X) = 1 - X(swap_X); 
end