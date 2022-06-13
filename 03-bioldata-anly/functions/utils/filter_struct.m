function s = filter_struct(s, f)
    % see: https://de.mathworks.com/matlabcentral/answers/1904-whats-a-fast-method-for-keeping-only-a-subset-of-fields-from-a-structure#answer_2842
    F = fieldnames(s); 
    D = struct2cell(s);
    M = cellfun(@(x) sum(contains(f, x)) > 0, fieldnames(s), 'uni', 1); 
    s = cell2struct(D(M,:), F(M));
end