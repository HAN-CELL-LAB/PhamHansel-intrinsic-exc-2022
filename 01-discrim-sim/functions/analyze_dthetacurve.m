function res = analyze_dthetacurve(dtheta_vec, anly_vec, sm_win)
if ~exist('sm_win', 'var'), sm_win = 1; end 
anly_vec = smooth(anly_vec, sm_win);
[pk, loc] = max(anly_vec); 

res.max = pk;
res.dtheta_at_max = dtheta_vec(loc); 
res.mean = mean(anly_vec);
res.mean_at_negdtheta = mean(anly_vec(dtheta_vec < 0));
res.mean_at_posdtheta = mean(anly_vec(dtheta_vec > 0));

end