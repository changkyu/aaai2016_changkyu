function show_results_cube( dir_result )

if( ~exist('dir_result','var') )
    env_setting;
end

experiments_desc = load_experiments_desc;
n_experiments = length(experiments_desc);

figure();
hold on;
for idx_ex=1:n_experiments

expr_desc = experiments_desc{idx_ex};

n_NVarr  = length(expr_desc.NVarr);
n_ETAarr = length(expr_desc.ETAarr);

n_params  = numel(expr_desc.params);
n_params2 = numel(expr_desc.params2);
n_run     = expr_desc.n_run;

%% Initialize
models_desc = load_models_desc;
n_models = numel(models_desc);
SAs = cell(n_models,1);
for idx_model=1:n_models
    model_desc = models_desc{idx_model};
    if sum(strcmp(model_desc.name, expr_desc.models)) == 0 
        continue; % skip if this experiment not including the current model
    end
        
    if( isequal(model_desc.type,'cppca') )
        SAs{idx_model} = zeros(n_params, n_params2);
    elseif( isequal(model_desc.type,'dppca') )
        SAs{idx_model} = cell(n_ETAarr, n_NVarr);
        for idx_ETA=1:n_ETAarr
            for idx_NV=1:n_NVarr
                NV  = expr_desc.NVarr(idx_NV);
                Networks = get_adj_graph(NV);
                n_Networks = length(Networks);                                
                SAs{idx_model}{idx_ETA,idx_NV} = cell(n_Networks,1);
                for idx_Network=1:n_Networks
                    SAs{idx_model}{idx_ETA,idx_NV}{idx_Network} ...
                     = zeros(n_params, n_params2);
                end
            end
        end
    end        
    
    tmp = load(fullfile( dir_result, [expr_desc.name '_' model_desc.name '.mat'] ));
    models_desc{idx_model}.models = tmp.models_desc{idx_model}.models;
end



%% Compute Subspace
for idx_ex_param=1:n_params
for idx_ex_param2=1:n_params2

    n_models = numel(models_desc);    
    for idx_model=1:n_models        
        model_desc = models_desc{idx_model};
        if sum(strcmp(model_desc.name, expr_desc.models)) == 0 
            continue; % skip if this experiment not including the current model
        end
            
        if isequal(model_desc.type, 'cppca')
            model = models_desc{idx_model}.models{idx_ex_param,idx_ex_param2};
            max_angle = subspace(model.W_GT, model.W);
            SAs{idx_model}(idx_ex_param,idx_ex_param2) = max_angle * 180 / pi;
        else
            for idx_ETA=1:n_ETAarr
                for idx_NV=1:n_NVarr
                    NV  = expr_desc.NVarr(idx_NV);
    
                    % Network topology
                    Networks = get_adj_graph(NV);
                    n_Networks = length(Networks);        
                    for idx_Network = 1:n_Networks
                        for idx_run = 1:n_run
    
                        model = models_desc{idx_model}.models{idx_ETA,idx_NV}{idx_Network}{idx_ex_param,idx_ex_param2};
                
                        max_angle = subspace(model.W_GT, model.W{1});    
                        for idx=2:NV
                            angle = subspace(model.W_GT, model.W{idx});
                            if max_angle < angle
                                max_angle = angle;
                            end
                        end
                        
                        SAs{idx_model}{idx_ETA,idx_NV}{idx_Network}(idx_ex_param,idx_ex_param2,idx_run)...
                         = max_angle * 180 / pi;
                     
                        end
                    end
                end
            end        
        end
    end    
end
end

%% Plot Boxes
hold on;
names_models = {'cppca', 'cppca', 'dppca', 'ours', 'ours (Tmax)'};
range_models = [2 4 5];
n_models = length(range_models);
colors = hsv(n_models);
x_step_size = ceil(n_models/3);
for idx=1:n_models
    idx_model = range_models(idx);
    
    position = x_step_size*idx_ex + ...
               (0.3*(-(n_models-1)/2)) + ...
               0.3*(idx-1);
    boxplot(SAs{idx_model}',...
        'labels',names_models{idx_model},...
        'colors',colors{idx,:},'symbol','g+','Position',position,'widths',0.2); 
    set(gca,'XTickLabel',{' '});
end
hold off;

% usually outliers are not too many...
% h=findobj(gca,'tag','Outliers');
% delete(h);

end
end