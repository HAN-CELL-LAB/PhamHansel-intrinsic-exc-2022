function new_vec = extendvec_fillmissing(vec_orig, len_new, fill_value)
if ~isvector(vec_orig)
    error('Only extend vectors');
end

if ~exist('fill_value', 'var'), fill_value = 0; end 

len_og = length(vec_orig);

if len_og >= len_new
    new_vec = vec_orig; 
    return;
end

sz = size(vec_orig); 
sz(find(sz == max(sz), 1)) = len_new;
new_vec = ones(sz) * fill_value; 
new_vec(1:len_og) = vec_orig;

end