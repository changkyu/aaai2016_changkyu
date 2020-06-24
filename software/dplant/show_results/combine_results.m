close all;
clearvars;

env_setting;

dir_result = '~/tmp2/results/';

models_desc = load_models_desc;
experiments_desc = load_experiments_desc;
range_experiments = [2 4];
range_models = [4:6 8:10 12:13];


for idx_ex=range_experiments
    expr_desc = experiments_desc{idx_ex};
    name_expr = expr_desc.name;
    fn_load_expr = expr_desc.fn_load;
    
    NVarr  = expr_desc.NVarr;
    ETAarr = expr_desc.ETAarr;
    n_NVarr  = length(expr_desc.NVarr);
    n_ETAarr = length(expr_desc.ETAarr);
    
    n_run = expr_desc.n_run;
    n_models = numel(models_desc);
    
    n_params  = numel(expr_desc.params );
    n_params2 = numel(expr_desc.params2);
    if n_params2==0
        n_params2 = 1;
    end
    
    results = cell(n_models,n_params,n_params2);
    for idx_model=range_models
        name_model = models_desc{idx_model}.name;        
        type_model = models_desc{idx_model}.type;
                
        if sum(strcmp(models_desc{idx_model}.name, expr_desc.models)) == 0 
            continue;
        end
    
        for idx_ex_param=1:n_params  
        params_expr = expr_desc.params;
        params2_expr = expr_desc.params2;

        for idx_ex_param2=1:n_params2

        if isempty(params2_expr)
            expr_params = params_expr{idx_ex_param};
        else
            expr_params = [ params_expr{idx_ex_param} 
                            params2_expr{idx_ex_param2} ];
        end
        
        if( strcmp(type_model,'cppca') )
            results{idx_model,idx_ex_param,idx_ex_param2}.iter = zeros(n_run,1);
            results{idx_model,idx_ex_param,idx_ex_param2}.time = zeros(n_run,1);
            %results{idx_model,idx_ex_param,idx_ex_param2}.err  = zeros(n_run,1);
            results{idx_model,idx_ex_param,idx_ex_param2}.err_rms  = zeros(n_run,1);
            results{idx_model,idx_ex_param,idx_ex_param2}.err_ssa  = zeros(n_run,1);
            results{idx_model,idx_ex_param,idx_ex_param2}.obj  =  cell(n_run,1);
            results{idx_model,idx_ex_param,idx_ex_param2}.ssa  =  cell(n_run,1);
            results{idx_model,idx_ex_param,idx_ex_param2}.rms  =  cell(n_run,1);            
            
        else
            results{idx_model,idx_ex_param,idx_ex_param2}.iter = zeros(n_ETAarr,n_NVarr,6,n_run);
            results{idx_model,idx_ex_param,idx_ex_param2}.time = zeros(n_ETAarr,n_NVarr,6,n_run);
            %results{idx_model,idx_ex_param,idx_ex_param2}.err  = zeros(n_ETAarr,n_NVarr,6,n_run);
            results{idx_model,idx_ex_param,idx_ex_param2}.err_rms  = zeros(n_ETAarr,n_NVarr,6,n_run);
            results{idx_model,idx_ex_param,idx_ex_param2}.err_ssa  = zeros(n_ETAarr,n_NVarr,6,n_run);
            results{idx_model,idx_ex_param,idx_ex_param2}.obj  =  cell(n_ETAarr,n_NVarr,6,n_run);            
            results{idx_model,idx_ex_param,idx_ex_param2}.ssa  =  cell(n_ETAarr,n_NVarr,6,n_run);            
            results{idx_model,idx_ex_param,idx_ex_param2}.rms  =  cell(n_ETAarr,n_NVarr,6,n_run);            
            results{idx_model,idx_ex_param,idx_ex_param2}.eta  =  cell(n_ETAarr,n_NVarr,6,n_run);            
        end

%         if( idx_ex==1 || idx_ex==2 )
%             [X, n_X, dim_X, dim_Z, GT] = fn_load_expr(expr_params{:});
%             fn_error = @calc_ppca_rms;
%         else
%             fn_error = @calc_ppca_max_ssa;
%         end

        for idx_run=1:n_run
            if( strcmp(type_model,'cppca') )
                save_mat = sprintf('%s_%s.mat',name_expr, name_model);
                filepath_save = fullfile(dir_result, save_mat);
                ld = load( filepath_save );
                model = ld.models;
                clear ld;

                if( strcmp(name_model,'cppca_affine') )
                    W_GT = model.W;
                else
                    if( strcmp(name_model,'cppca_em'))
                        results{idx_model,idx_ex_param,idx_ex_param2}.iter(idx_run) = model.eITER;
                        results{idx_model,idx_ex_param,idx_ex_param2}.time(idx_run) = model.eTIME;
                    end
%                     if( isequal(fn_error,@calc_ppca_rms) )
%                         err = fn_error(X, model.W, model.EZ, model.MU);
%                     else
%                         err = fn_error(model.W_GT, model.W);
%                     end
                    [X, n_X, dim_X, dim_Z, GT] = fn_load_expr(expr_params{:});
                    err_rms = calc_ppca_rms(X, model.W, model.EZ, model.MU);
                    err_ssa = calc_ppca_max_ssa(model.W_GT, model.W);
                    
                    results{idx_model,idx_ex_param,idx_ex_param2}.err_rms(idx_run) = err_rms;
                    results{idx_model,idx_ex_param,idx_ex_param2}.err_ssa(idx_run) = err_ssa;
                    results{idx_model,idx_ex_param,idx_ex_param2}.obj{idx_run} = model.objArray;
                    results{idx_model,idx_ex_param,idx_ex_param2}.ssa{idx_run} = model.ssaArray;                    
                    results{idx_model,idx_ex_param,idx_ex_param2}.rms{idx_run} = model.rmsArray;                                    
                end
                clear model;
            else
                for idx_ETA=1:n_ETAarr
                    ETA = ETAarr(idx_ETA);
                    for idx_NV=1:n_NVarr
                        NV = NVarr(idx_NV);
                        Networks = get_adj_graph_ex(NV);
                        n_Networks = length(Networks);
                        for idx_Network=1:n_Networks
                            Network = Networks{idx_Network};

                            save_mat = sprintf('%s_%s_ETA%d_NV%d_%s.mat',name_expr,name_model,ETA,NV,Network.name);
                            filepath_save = fullfile(dir_result, save_mat);
                            ld = load( filepath_save );
                            model = ld.models{idx_ex_param,idx_ex_param2,idx_run};
                            if( strcmp(name_model,'dppca_ant') || strcmp(name_model,'dppca_ant_Tmax') )
                                if( size(model.ETA_ij_history,1) > model.eITER )
                                    models = ld.models;

                                    for idx_p1=1:n_params
                                    for idx_p2=1:n_params2
                                    for idx_r=1:n_run

                                    models{idx_p1,idx_p2,idx_r}.ETA_ij_history = models{idx_p1,idx_p2,idx_r}.ETA_ij_history(1:models{idx_p1,idx_p2,idx_r}.eITER,:,:);
                                    models{idx_p1,idx_p2,idx_r}.tau_ij_history = models{idx_p1,idx_p2,idx_r}.tau_ij_history(1:models{idx_p1,idx_p2,idx_r}.eITER,:,:);
                                    models{idx_p1,idx_p2,idx_r}.T_ij_history   = models{idx_p1,idx_p2,idx_r}.T_ij_history(1:models{idx_p1,idx_p2,idx_r}.eITER,:,:);
                                    models{idx_p1,idx_p2,idx_r}.comm_histroy   = models{idx_p1,idx_p2,idx_r}.comm_histroy(1:models{idx_p1,idx_p2,idx_r}.eITER,:,:);

                                    end
                                    end
                                    end

                                    save(filepath_save,'models','-v7.3');
                                    disp(filepath_save);

                                    model = ld.models{idx_ex_param,idx_ex_param2,idx_run};
                                    clear models;
                                end
                            end
                            clear ld;

                            results{idx_model,idx_ex_param,idx_ex_param2}.iter(idx_ETA,idx_NV,idx_Network,idx_run) = model.eITER;
                            results{idx_model,idx_ex_param,idx_ex_param2}.time(idx_ETA,idx_NV,idx_Network,idx_run) = model.eTIME;
%                             if( isequal(fn_error,@calc_ppca_rms) )
%                                 tmp = zeros(NV,1);
%                                 for i=1:NV
%                                     tmp(i) = fn_error(X, model.W, model.EZ, model.MU);
%                                 end
%                                 
%                                 err = max(tmp);
%                             else
%                                 if( ~isempty(model.W_GT) )
%                                     W_GT = model.W_GT;
%                                     tmp = zeros(NV,1);
%                                     for i=1:NV
%                                         tmp(i) = fn_error(W_GT, model.W{i});
%                                     end
%                                     
%                                     err = max(tmp);
%                                 else
%                                     err = model.ssaArray(end);
%                                 end                                
%                             end
                            err_ssa = model.ssaArray(end);
                            err_rms = model.rmsArray(end);

                            results{idx_model,idx_ex_param,idx_ex_param2}.err_rms(idx_ETA,idx_NV,idx_Network,idx_run) = err_rms;
                            results{idx_model,idx_ex_param,idx_ex_param2}.err_ssa(idx_ETA,idx_NV,idx_Network,idx_run) = err_ssa;
                            results{idx_model,idx_ex_param,idx_ex_param2}.obj{idx_ETA,idx_NV,idx_Network,idx_run} = model.objArray;
                            results{idx_model,idx_ex_param,idx_ex_param2}.ssa{idx_ETA,idx_NV,idx_Network,idx_run} = model.ssaArray;
                            results{idx_model,idx_ex_param,idx_ex_param2}.rms{idx_ETA,idx_NV,idx_Network,idx_run} = model.rmsArray;                            
                            if( strcmp(name_model,'dppca_ant') || strcmp(name_model,'dppca_ant_Tmax'))
                                results{idx_model,idx_ex_param,idx_ex_param2}.eta{idx_ETA,idx_NV,idx_Network,idx_run} = model.ETA_ij_history;
                            else
                                results{idx_model,idx_ex_param,idx_ex_param2}.eta{idx_ETA,idx_NV,idx_Network,idx_run} = ETA*ones(model.eITER,NV,NV);
                            end
                            
                            clear model;
                        end
                    end
                end
            end
        end
        end
        end
    end
    
    if( ~exist( fullfile(dir_result, 'summary'), 'dir' ) )
        mkdir(fullfile(dir_result, 'summary'));
    end
    save(fullfile(dir_result, 'summary', sprintf('%s_summary.mat',name_expr)), 'results','-v7.3');
    clear results;    
end

