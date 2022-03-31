function results = get_specperf(opts)

demo_spec = false; 
if isfield(opts, 'demo') 
    demo_spec = opts.demo; 
end

if isfield(opts, 'num_Wrand')
    num_Wrand = opts.num_Wrand; 
    opts = rmfield(opts, 'num_Wrand'); 
    tmp_res = arrayfun(@(~) get_specperf(opts), 1:num_Wrand, 'uni', 1); 
    
    if demo_spec
        results = tmp_res; 
        return;
    end
    
    res_fields = fieldnames(tmp_res); 
    results = struct; 
    for i = 1:length(res_fields)
        results.(res_fields{i}) = mean(horzcat(tmp_res.(res_fields{i})), 2);
    end
    return; 
end

N_Y = opts.N_Y;
N_X = opts.N_X;
alpha_W = opts.alpha_W;
num_dthres = opts.num_dthres;
thres_baseline = opts.thres_baseline;
num_inputs = opts.num_inputs;
sigma_noise = opts.sigma_noise;
rng_dthres = [-0.5, 0.5]; 
if isfield(opts, 'rng_dthres')
    rng_dthres = opts.rng_dthres; 
end

percent_complete_input = 1; 
if isfield(opts, 'percent_complete_input')
    percent_complete_input = opts.percent_complete_input; 
end

num_overlap_per_group = 0; 
if isfield(opts, 'num_overlap_per_group') 
    num_overlap_per_group = opts.num_overlap_per_group;
end

N_XselperY = round(N_X/N_Y);
W_YX = zeros(N_Y, N_X);
dthres_vec = linspace(rng_dthres(1),rng_dthres(2),num_dthres);

selective_indices = arrayfun(@(i) (i-1)*N_XselperY+1:i*N_XselperY, 1:N_Y, 'uni', 0);
if num_overlap_per_group > 0
    overlap_inds = 1:num_overlap_per_group;
    N_XselperY_nonoverlap = N_XselperY - num_overlap_per_group; 
    selective_indices = arrayfun(@(i) [overlap_inds, num_overlap_per_group + ((i-1)*N_XselperY_nonoverlap+1:i*N_XselperY_nonoverlap)], 1:N_Y, 'uni', 0);
end

num_present_input = round(N_XselperY * percent_complete_input); 

representative_inputs = cell(N_Y,1);
for i = 1:N_Y
    W_YiX = rand(N_X,1);
    sel_ind = selective_indices{i};
    nonsel_ind = true(N_X,1);
    nonsel_ind(sel_ind) = false;
    
    W_YiX(sel_ind) = rand(N_XselperY,1);
    W_YiX(sel_ind) = alpha_W * W_YiX(sel_ind)/sum(W_YiX(sel_ind));
    
    W_YiX(nonsel_ind) = (1-alpha_W) *  W_YiX(nonsel_ind) / sum(W_YiX(nonsel_ind));
    
    W_YX(i,:) = W_YiX;
    
    X_of_Yi = zeros(N_X,1);
    X_of_Yi(sel_ind) = 1;
    representative_inputs{i} = X_of_Yi;
end

if demo_spec
    results = struct;
    labels = repmat(1:N_Y, [1,num_inputs]); 
    raw_inputs = representative_inputs(labels);
    results.labels = labels; 
    results.raw_inputs = raw_inputs;
    results.W_YX = W_YX;
    
    masked_inputs = cellfun(@(x) mask_input(x, num_present_input), raw_inputs, 'uni', 0);
    background_noise = cellfun(@(x) sigma_noise * randn(size(x)), masked_inputs, 'uni', 0);
    actual_inputs = cellfun(@(x,y) x+y, masked_inputs, background_noise, 'uni', 0);
    
    results.masked_inputs = masked_inputs;
    results.background_noise = background_noise;
    results.actual_inputs = actual_inputs;

    return; 
end

labels = randi(N_Y, [1,num_inputs]);
Xs = representative_inputs(labels);
Xs = cellfun(@(x) mask_input(x, num_present_input), Xs, 'uni', 0); 
Xs = horzcat(Xs{:});
Xs = Xs + sigma_noise * randn(size(Xs));

preY = W_YX * Xs;

discrim_acc_vec = zeros(num_dthres,1);
recog_acc_vec = zeros(num_dthres,1); 

recog_TPR_vec = zeros(num_dthres,1);
recog_FPR_vec = zeros(num_dthres,1);
recog_TNR_vec = zeros(num_dthres,1);
recog_FNR_vec = zeros(num_dthres,1);

thresY_basevec = thres_baseline * ones(N_Y,1);
label_onehot_mat = labels  == (1:N_Y)'; 

for i = 1:num_dthres
    tmp_discrim = zeros(N_Y,1);
    tmp_recog = zeros(N_Y,1); 
    
    tmp_recog_TPR = zeros(N_Y,1);
    tmp_recog_FPR = zeros(N_Y,1);
    tmp_recog_TNR = zeros(N_Y,1);
    tmp_recog_FNR = zeros(N_Y,1);
    
    for k = 1:N_Y
        thres_Y = thresY_basevec;
        thres_Y(k) = thres_Y(k) + dthres_vec(i);
        Y = preY > thres_Y;        
        
        % tmp_acc(k) = mean(arrayfun(@(x) sum(Y(:,x)) == 1 && Y(labels(x),x), 1:num_inputs, 'uni', 1));
        % % below faster, test by: 
        %   ninp = 2000; nlbl = 10; yy = logical(randi(2,nlbl,ninp) - 1); lbl = randi(nlbl,1,ninp);
        %   tic; abc_1 = arrayfun(@(x) sum(yy(:,x)) == 1 && yy(lbl(x),x), 1:ninp, 'uni', 1); toc
        %   tic; lblmat = lbl == (1:nlbl)'; abc_2 = sum(yy,1) == 1 & sum(abs(yy-lblmat),1) == 0;  toc
        %   isequal(abc_1, abc_2)        
        tmp_discrim(k) = mean(sum(Y,1) == 1 & sum(abs(Y - label_onehot_mat),1) == 0);
        
        tmp_recog(k) = sum(Y(k,:) == 1 & labels == k) / sum(labels == k);
        
        real_labels_of_Yk = label_onehot_mat(k,:)';
        pred_labels_of_Yk = Y(k,:)'; 
        tmp_recogperfrate = classif_stats(real_labels_of_Yk, pred_labels_of_Yk);
        
        tmp_recog_TPR(k) = tmp_recogperfrate.TPR;
        tmp_recog_FPR(k) = tmp_recogperfrate.FPR;
        tmp_recog_TNR(k) = tmp_recogperfrate.TNR;
        tmp_recog_FNR(k) = tmp_recogperfrate.FNR;
        
    end
    
    discrim_acc_vec(i) = mean(tmp_discrim);
    recog_acc_vec(i) = mean(tmp_recog); 
    
    recog_TPR_vec(i) = mean(tmp_recog_TPR);
    recog_FPR_vec(i) = mean(tmp_recog_FPR);
    recog_TNR_vec(i) = mean(tmp_recog_TNR);
    recog_FNR_vec(i) = mean(tmp_recog_FNR);
    
end

results = struct; 
results.discrimination_accuracy = discrim_acc_vec;
results.recognition_accuracy = recog_acc_vec;

results.recognition_TPR = recog_TPR_vec;
results.recognition_FPR = recog_FPR_vec;
results.recognition_TNR = recog_TNR_vec;
results.recognition_FNR = recog_FNR_vec;

end

function X = mask_input(X, num_present)
ind_on = find(X); 
X(ind_on(randperm(length(ind_on), length(ind_on) - num_present))) = 0;  
end