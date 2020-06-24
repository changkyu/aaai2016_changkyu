% Experiment for Distributed Prob. Learning with Adaptive Network Topology
%--------------------------------------------------------------------------
%
% Implemented/Modified from [1]
%  by     Changkyu Song (changkyu.song@rutgers.edu)
%  on     2015.03.19 (last modified on 2015.03.19)
%
% References
%  [1] S. Yoon and V. Pavlovic. Distributed Probabilistic Learning
%      for Camera Networks with Missing Data. In NIPS, 2012.
%
%--------------------------------------------------------------------------

%% Experiments Description
experiments_desc = load_experiments_desc;
models_desc = load_models_desc;

if( ~exist('range_experiments', 'var') )
    range_experiments = 1:numel(experiments_desc);
end
if( ~exist('range_models', 'var') )
    range_models = 1:numel(models_desc);
end
b_set_Networks = exist('range_Networks', 'var');

for idx_ex=range_experiments%1:length(experiments_desc)

expr_desc          = experiments_desc{idx_ex};
name_expr          = expr_desc.name;
fmt_name_expr      = expr_desc.fmt_name;
fmt_result_expr    = expr_desc.fmt_result;
fn_load_expr       = expr_desc.fn_load;
params_expr        = expr_desc.params;
params2_expr       = expr_desc.params2;
b_zeromean_expr    = expr_desc.ZeroMean;
init_params_c_expr = expr_desc.init_params_c;
init_params_d_expr = expr_desc.init_params_d;
objfreq_expr       = expr_desc.objfreq; % freq for printing iter
objprec_expr       = expr_desc.objprec; % prec for convergence
type_expr          = expr_desc.type;
n_run_expr         = expr_desc.n_run;

n_params = numel(params_expr);
n_params2 = numel(params2_expr);
if n_params2==0
    n_params2 = 1;
end

% NV: the number of vertices (nodes) in the network
n_NVarr = length(expr_desc.NVarr);

% ETA: Learning rate
n_ETAarr = length(expr_desc.ETAarr);

%% Check dataset availability
if strfind(name_expr, 'cube')
    dir_data = fullfile(dir_dataset,'cube');
    if ~exist(dir_data, 'dir')
        synthetic_cube_data_gen_ex(dir_dataset);
    end    
elseif strfind(name_expr, 'caltech')
    dir_data = fullfile(dir_dataset,'caltech');
    if ~exist(dir_data, 'dir')
        error('Please put Caltech Turntable Dataset into data the directory.');
    end    
elseif strfind(name_expr, 'hopkins')
    dir_data = fullfile(dir_dataset,'hopkins');
    if ~exist(dir_data, 'dir')
        error('Please put Hopkins 155 Dataset into data the directory.');
    end    
end

%% Experiments - Centralized
disp(['============ ' name_expr ' ============']);

for idx_ex_param=1:n_params        
    for idx_ex_param2=1:n_params2
        if isempty(params2_expr)
            expr_params = params_expr{idx_ex_param};
        else
            expr_params = [ params_expr{idx_ex_param} 
                            params2_expr{idx_ex_param2} ];
        end

    	fprintf(['Expr Detail: ' sprintf(fmt_name_expr, expr_params{:}) '\n']);

        [X, n_X, dim_X, dim_Z, GT] = feval(fn_load_expr, expr_params{:});
        if b_zeromean_expr
            X_zm = X - repmat(mean(X,1), [n_X,1]);
        end

        init_model = get_init_value_ex(X, dim_Z, init_params_c_expr{:});

        for idx_model=range_models
                
            if isequal(models_desc{idx_model}.type, 'dppca')
                continue; % skip distributed models
            end
            if sum(strcmp(models_desc{idx_model}.name, expr_desc.models)) == 0 
                continue; % skip if this experiment not including the current model
            end

            desc_model       = models_desc{idx_model}.desc;
            fn_model         = models_desc{idx_model}.fn;
            name_model       = models_desc{idx_model}.name;

            disp(['============ ' name_model ' ============']);
            
            models    = cell(n_params, n_params2, n_run_expr);
            parfor idx_run=1:n_run_expr            
                    models{idx_ex_param, idx_ex_param2,idx_run} ...
                     = parfor_cppca( desc_model, fn_model, X, dim_Z, objprec_expr, objfreq_expr, init_model);
                    models{idx_ex_param, idx_ex_param2,idx_run}.W_GT = GT; % ground truth
            end
            save( fullfile(dir_result, [ expr_desc.name '_' model_desc.name '.mat']), 'model_desc', 'models', '-v7.3');
        end
    end
end

%% Experiments - Distributed
if( n_run_expr==1 )                
    
    META = cell(n_ETAarr*n_NVarr*6*n_params*n_params2,1);
    idx_META = 0;
    for idx_ETA = 1:n_ETAarr
        ETA = expr_desc.ETAarr(idx_ETA);
        for idx_NV = 1:n_NVarr
            NV = expr_desc.NVarr(idx_NV);
            Networks = get_adj_graph_ex(NV);
            n_Networks = length(Networks);        
            if( b_set_Networks==false )
                range_Networks = 1:n_Networks;
            end
            for idx_Network = range_Networks
                for idx_ex_param=1:n_params
                for idx_ex_param2=1:n_params2
                for idx_run=1:n_run_expr

                    if isempty(params2_expr)
                        expr_params = params_expr{idx_ex_param};
                    else
                        expr_params = [ params_expr{idx_ex_param} 
                                        params2_expr{idx_ex_param2} ];
                    end

                    idx_META = idx_META + 1;
                    META{idx_META} = {{ idx_ETA,     ...
                                        ETA,         ...
                                        idx_NV,      ...
                                        NV,          ...
                                        idx_Network, ...
                                        Networks{idx_Network}, ...
                                        idx_ex_param,idx_ex_param2, ...
                                        expr_params, idx_run }};
                end
                end
                end
            end                
        end
    end
    n_META = idx_META;        
    models_meta = cell(n_ETAarr, n_NVarr);
    for idx_ETA = 1:n_ETAarr
        ETA = expr_desc.ETAarr(idx_ETA);
        for idx_NV = 1:n_NVarr
            NV = expr_desc.NVarr(idx_NV);
            Networks = get_adj_graph_ex(NV);
            n_Networks = length(Networks);        
            models_meta{idx_ETA,idx_NV} = cell(n_Networks,1);
            if( b_set_Networks==false )
                range_Networks = 1:n_Networks;
            end
            for idx_Network = range_Networks
                models_meta{idx_ETA,idx_NV}{idx_Network} ...
                 = cell(n_params,n_params2,n_run_expr);
            end
        end
    end

    parfor idx_meta=1:n_META
        meta = META{idx_meta}{:};
        idx_ETA       = meta{1};
        ETA           = meta{2};
        idx_NV        = meta{3};
        NV            = meta{4};
        idx_Network   = meta{5};
        name_Network  = meta{6}.name;
        adj_Network   = meta{6}.adj;
        idx_ex_param  = meta{7};
        idx_ex_param2 = meta{8};
        expr_params   = meta{9};
        idx_run       = meta{10};            

        [X, n_X, dim_X, dim_Z, GT] = feval(fn_load_expr,expr_params{:});
        if b_zeromean_expr
            X = X - repmat(mean(X,1), [n_X,1]);
        end

        fprintf('NV: %d, dim_X: %d\n', NV, dim_X);
        if ( isequal(type_expr,'sfm') )
            Vp = get_sample_assign(NV, dim_X, 'IsSfM', true);
        else    
            Vp = get_sample_assign(NV, n_X);
        end
        init_model = get_init_value_ex( X, dim_Z,  ...
                     [init_params_d_expr,          ...
                     {'NumberOfNodes',   NV, 'SampleAssignVec', Vp}]);

        disp(['idx_meta: ', num2str(idx_meta)]);
        
        for idx_model=range_models
            
            if isequal(model_desc.type, 'cppca')
                continue; % skip centralized models
            end
            if sum(strcmp(model_desc.name, expr_desc.models)) == 0
                continue; % skip if this experiment not including the current model
            end
        
            disp(['============ ' models_desc{idx_model}.name ' ============']);

            model ...
             = parfor_dppca(models_desc{idx_model},      ...
                            X,dim_Z,GT,b_zeromean_expr, ...
                            ETA,Vp,adj_Network,               ...
                            objprec_expr,objfreq_expr,init_model);                      
            model.W_GT = GT; % ground truth
            disp(['Iter: ' num2str(model.eITER)]);
            save_parfor_model( fullfile(dir_result, sprintf('%s_%s_ETA%d_NV%d_%s_idxp%d_idxp2%d_idxr%d.mat', ...
                               name_expr, models_desc{idx_model}, ETA, NV, name_Network, idx_ex_param, idx_ex_param2, idx_run)), model);        
        end
    end

    for idx_model=range_models    
        for idx_ETA = 1:n_ETAarr
            ETA = expr_desc.ETAarr(idx_ETA);
            for idx_NV = 1:n_NVarr
                NV = expr_desc.NVarr(idx_NV);
                Networks = get_adj_graph_ex(NV);
                n_Networks = length(Networks);        
                models_meta{idx_ETA,idx_NV} = cell(n_Networks,1);
                if( b_set_Networks==false )
                    range_Networks = 1:n_Networks
                end
                for idx_Network = range_Networks
                    name_Network = Networks{idx_Network}.name;
                    models = cell(n_params, n_params2, n_run_expr);
                    for idx_ex_param=1:n_params
                    for idx_ex_param2=1:n_params2
                    for idx_run=1:n_run_expr
                        filepath_load = fullfile(dir_result, sprintf('%s_%s_ETA%d_NV%d_%s_idxp%d_idxp2%d_idxr%d.mat', ...
                               name_expr, name_model, ETA, NV, name_Network, idx_ex_param, idx_ex_param2, idx_run));
                        ld = load( filepath_load );
                        models{idx_ex_param,idx_ex_param2,idx_run} = ld.model;
                        clear ld;
                        system(['rm ' filepath_load]);
                    end
                    end
                    end
                    save( fullfile(dir_result, sprintf('%s_%s_ETA%d_NV%d_%s.mat', expr_desc.name, models_desc{idx_model}.name, ETA, NV, name_Network)), 'model_desc', 'models', '-v7.3');
                    clear models;
                end
            end
        end    
    end
else

    for idx_ETA = 1:n_ETAarr
        ETA = expr_desc.ETAarr(idx_ETA);
        fprintf(['ETA: ' num2str(ETA) '\n']);

        for idx_NV = 1:n_NVarr
            % Node assignment to each sample
            NV = expr_desc.NVarr(idx_NV);
            fprintf(['# of nodes: ' num2str(NV) '\n']);

            % Network topology
            Networks = get_adj_graph_ex(NV);
            n_Networks = length(Networks);        
            if( b_set_Networks==false )
                range_Networks = 1:n_Networks
            end
            for idx_Network = range_Networks
                name_Network = Networks{idx_Network}.name;
                adj_Network  = Networks{idx_Network}.adj;
                fprintf(['Network: ' name_Network '\n']);

                for idx_ex_param=1:n_params
                    for idx_ex_param2=1:n_params2
                        if isempty(params2_expr)
                            expr_params = params_expr{idx_ex_param};
                        else
                            expr_params = [ params_expr{idx_ex_param} 
                                            params2_expr{idx_ex_param2} ];
                        end

                        fprintf(['Expr Detail: ' sprintf(fmt_name_expr, expr_params{:}) '\n']);

                        [X, n_X, dim_X, dim_Z, GT] = fn_load_expr(expr_params{:});
                        if b_zeromean_expr
                            X = X - repmat(mean(X,1), [n_X,1]);
                        end
                        
                        fprintf('NV: %d, dim_X: %d\n', NV, dim_X);
                        if ( isequal(type_expr,'sfm') )
                            %Vp = get_sample_assign(NV, dim_X, 'IsSfM', true);
                            Vp = get_sample_assign(NV, dim_X);
                        else    
                            Vp = get_sample_assign(NV, n_X);
                        end
                        
                        parfor idx_run=1:n_run_expr
                        
                            init_model = get_init_value_ex( X, dim_Z,  ...
                                         [init_params_d_expr,          ...
                                         {'NumberOfNodes',   NV, 'SampleAssignVec', Vp}]);

                            for idx_model=range_models

                                if isequal(models_desc{idx_model}.type, 'cppca')
                                    continue; % skip centralized models
                                end
                                if sum(strcmp(models_desc{idx_model}.name, expr_desc.models)) == 0
                                    continue; % skip if this experiment not including the current model
                                end

                                model ...
                                 = parfor_dppca(models_desc{idx_model},     ...
                                                X,dim_Z,GT,b_zeromean_expr, ...
                                                ETA,Vp,adj_Network,         ...
                                                objprec_expr,objfreq_expr,init_model);                      
                                model.W_GT = GT; % ground truth
                                save_parfor_model( fullfile(dir_result, sprintf('%s_%s_ETA%d_NV%d_%s_idxp%d_idxp2%d_idxr%d.mat', ...
                                                   name_expr, models_desc{idx_model}.name, ETA, NV, name_Network, idx_ex_param, idx_ex_param2, idx_run)), model);
                               
                                fprintf('(Iter: %d,\tssa(o): %.2f\t) %s\n',model.eITER,model.ssaArray(model.eITER)/pi*180, models_desc{idx_model}.name);

                            end
                        end
                    end                
                end

            end
        end
    end
    
    for idx_model=range_models    
        model_desc = models_desc{idx_model};
        for idx_ETA = 1:n_ETAarr
            ETA = expr_desc.ETAarr(idx_ETA);
            for idx_NV = 1:n_NVarr
                NV = expr_desc.NVarr(idx_NV);
                Networks = get_adj_graph_ex(NV);
                n_Networks = length(Networks);        
                models_meta{idx_ETA,idx_NV} = cell(n_Networks,1);
                if( b_set_Networks==false )
                    range_Networks = 1:n_Networks;
                end
                for idx_Network = range_Networks
                    name_Network = Networks{idx_Network}.name;
                    models = cell(n_params, n_params2, n_run_expr);
                    for idx_ex_param=1:n_params
                    for idx_ex_param2=1:n_params2
                    for idx_run=1:n_run_expr
                        filepath_load = fullfile(dir_result, sprintf('%s_%s_ETA%d_NV%d_%s_idxp%d_idxp2%d_idxr%d.mat', ...
                               name_expr, model_desc.name, ETA, NV, name_Network, idx_ex_param, idx_ex_param2, idx_run));
                        ld = load( filepath_load );
                        models{idx_ex_param,idx_ex_param2,idx_run} = ld.model;
                        clear ld;
                        system(['rm ' filepath_load]);
                    end
                    end
                    end
                    save( fullfile(dir_result, sprintf('%s_%s_ETA%d_NV%d_%s.mat', expr_desc.name, models_desc{idx_model}.name, ETA, NV, name_Network)), 'model_desc', 'models', '-v7.3');
                    clear models;
                end
            end
        end    
    end
end

end
