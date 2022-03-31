function res = classif_stats(real_labels, pred_labels)
% the dim of labels must be: [num_cases x num_conditions]

real_labels = real_labels > 0;
pred_labels = pred_labels > 0;

% simple stats 
num_real_P = sum(real_labels,1);
num_real_N = sum(~real_labels,1);

TP = sum(real_labels  &  pred_labels,1);  % true positives
TN = sum(~real_labels &  ~pred_labels,1); % true negatives

FP = sum(~real_labels &  pred_labels,1);  % false positives
FN = sum(real_labels  &  ~pred_labels,1); % false negatives

res.TP = TP;
res.TN = TN; 
res.FP = FP;
res.FN = FN; 

% rates 
TPR = TP ./ num_real_P; % true positive rate
FNR = FN ./ num_real_P; % false negative rate

FPR = FP ./ num_real_N; % false positive rate
TNR = TN ./ num_real_N; % true negative rate 

res.TPR = TPR;
res.FNR = FNR;
res.FPR = FPR; 
res.TNR = TNR;

% other stats/aliases
res.accuracy = (TP + TN) ./ (num_real_P + num_real_N); 

res.sensitity = TPR; 
res.recall = TPR; 
res.hitrate = TPR; 

res.specificity = TNR; 
res.selectivity = TNR; 

res.precision = TP ./ (TP + FP);
end

