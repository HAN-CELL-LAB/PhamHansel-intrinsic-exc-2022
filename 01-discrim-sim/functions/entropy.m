function ent = entropy(p)
    p = p / sum(p); 
    ent = -p .* log2(p); 
    ent = sum(ent(~isinf(ent) & ~isnan(ent))); 
end