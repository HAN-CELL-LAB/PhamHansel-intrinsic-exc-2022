classdef classifyObjA < handle
    
    properties
        % passed into constructor
        N_side;
        N_Y;
        
        size_input_dataset;
        prop_objA_in_dataset;
        
        num_thres;
        
        name_of_objA;
        args_of_objA = {};
        
        name_of_notA;
        args_of_notA = {};
                
        genWmethod;
        
        genWmain_name;
        genWmainfun;
        
        genWnoise_name;
        genWnoisefun;
        alpha_W;
        
        % defined based on constructor
        N_X;
        thres_vec;
        X_A;
        W_YX;
        num_unqactiveX; 
        
        % example dataset (instead of saving all) 
        example_dataset = struct('objA_full', '', 'objA_input', '', 'notA_full', '', 'notA_input', ''); 
        
        % output
        stats = struct('fields', '', 'all', '', 'basedon_numactiveX', ''); 
    end
    
    methods
        function obj = classifyObjA(props)
            prop_names = fieldnames(props);
            for ind_prop = 1:length(prop_names)
                prop_name = prop_names{ind_prop};
                obj.(prop_name) = props.(prop_name);
            end
            
            obj.N_X = obj.N_side^2;
            obj.thres_vec = linspace(0,1,obj.num_thres);
            obj.X_A = to_col_vec(create_squareobj(obj.N_side, obj.name_of_objA, obj.args_of_objA{:}));            
            obj.num_unqactiveX = 0:obj.N_X;            
            
        end
        
        function run(obj)
            generate_weights(obj); % define W_YX
            
            evaluate_outputstats(obj); % stats 
        end
        
        function evaluate_outputstats(obj)
            % Create dataset
            [Xs, real_labels] = obj.create_dataset; 
            num_activeX = sum(Xs, 1);
            
            % Output based on thresholds
            Y = obj.W_YX * Xs;
            Y = Y' > obj.thres_vec;
            
            % Stats 
            pred_labels = Y;
                        
            obj.stats.all = classif_stats(real_labels, pred_labels);
            
            stats_fields = fieldnames(obj.stats.all);
            obj.stats.fields = stats_fields; 
            
            % States based on number of active elements in X
            stats_basedon_numactveX = arrayfun(@(k) ...
                classif_stats(real_labels(num_activeX == k), pred_labels(num_activeX == k, :)), ...
                obj.num_unqactiveX, 'uni', 0);
            stats_basedon_numactveX = vertcat(stats_basedon_numactveX{:});
            stats_basedon_numactveX = cellfun(@(n) vertcat(stats_basedon_numactveX.(n)), stats_fields, 'uni', 0);
            
            obj.stats.basedon_numactiveX = cell2struct(stats_basedon_numactveX, stats_fields);
            
        end
        

        function [Xs_masked, labels] = create_dataset(obj)
            
            num_objA_in_dataset = round(obj.size_input_dataset * obj.prop_objA_in_dataset);
            num_notA_in_dataset = obj.size_input_dataset - num_objA_in_dataset;
            
            labels = zeros(obj.size_input_dataset,1);
            labels(1:num_objA_in_dataset) = 1;
            
            Xs = zeros(obj.N_X, obj.size_input_dataset);
            
            % creating objA's
            Xs(:,labels == 1) = repmat(obj.X_A, [1, num_objA_in_dataset]);
            
            % creating notA's
            if ~strcmpi(obj.name_of_notA, 'NOISE')
                X_notA = to_col_vec(create_squareobj(obj.N_side, obj.name_of_notA, obj.args_of_notA{:}));
                Xs(:, labels == 0) = repmat(X_notA, [1, num_notA_in_dataset]);
            else
                % arrayfun is much slower, so only used with "NOISE" generation
                Xs_notA = arrayfun(@(~) to_col_vec(create_squareobj(obj.N_side, 'NOISE', obj.args_of_notA{:})), ...
                    1:num_notA_in_dataset, 'uni', 0);
                Xs(:, labels == 0) = horzcat(Xs_notA{:});
            end
            
            % shuffle labels and dataset
            shuffled_inds = randperm(obj.size_input_dataset);
            Xs = Xs(:,shuffled_inds);
            labels = labels(shuffled_inds);
            
            % masking (covering) inputs 
            mask_Xs = rand(obj.N_X, obj.size_input_dataset);
            p_notcovered = (randi(obj.N_X+1,[1,obj.size_input_dataset]) - 1)/obj.N_X;
            mask_Xs = mask_Xs < p_notcovered;
            Xs_masked = Xs .* mask_Xs;
            
            % save examples 
            sz_sq = [1,1]*obj.N_side; 
            ind_example_objA = find(labels==1,1);
            ind_example_notA = find(labels==0,1);
            
            obj.example_dataset.objA_full = reshape(obj.X_A, sz_sq);
            obj.example_dataset.objA_input = reshape(Xs_masked(:,ind_example_objA), sz_sq);
            
            obj.example_dataset.notA_full = reshape(Xs(:,ind_example_notA), sz_sq);
            obj.example_dataset.notA_input = reshape(Xs_masked(:,ind_example_notA), sz_sq);
        
        end
        
        function generate_weights(obj)
            
            switch upper(obj.genWmethod)
                case 'ALL'
                    W = obj.genWmainfun(obj.N_Y,obj.N_X);
                    W = W/sum(W);
                case 'SEPARATE'
                    ind_of_objA = obj.X_A == 1; 
                    ind_oj_notA = ~ind_of_objA;
                    
                    n_X_objA = sum(ind_of_objA);
                    n_X_notA = obj.N_X - n_X_objA;
                    
                    W = zeros(obj.N_Y,obj.N_X);
                    
                    W(:,ind_of_objA) = obj.genWmainfun(obj.N_Y,n_X_objA);
                    W(:,ind_of_objA) = obj.alpha_W * W(:,ind_of_objA)/sum(W(:,ind_of_objA),'all');
                    
                    W(:,ind_oj_notA) = obj.genWnoisefun(obj.N_Y,n_X_notA);
                    W(:,ind_oj_notA) = (1-obj.alpha_W) * W(:,ind_oj_notA)/sum(W(:,ind_oj_notA),'all');
                otherwise
                    error('"%s" is not an allowed weight generation method.', obj.genWmethod);
            end
            
            obj.W_YX = W;
        end
    end
end