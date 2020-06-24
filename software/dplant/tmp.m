experiments_desc = load_experiments_desc;
n_experiments = numel(experiments_desc);

for idx_expr=1:n_experiments
    expr_desc = experiments_desc{idx_expr};
    names_models = expr_desc.models;
    n_models = numel(names_models);
    
    ETAarr = expr_desc.ETAarr;
    n_ETAarr  = numel(ETAarr);
    
    NVarr = expr_desc.NVarr;
    n_NVarr  =numel(NVarr);
    for idx_m=1:n_models
        name_model = names_models{idx_m};
        switch name_model
            case 'cppca_em'
                idx_model = 2;
            case 'cppca_affine'
                idx_model = 3;
            case 'dppca'
                idx_model = 4;
            case 'dppca_ant'
                idx_model = 5;
            case 'dppca_ant_Tmax'
                idx_model = 6;
        end

        name_mat = sprintf('%s_%s.mat',expr_desc.name,names_models{idx_m});
        disp(name_mat);
        ld = load(fullfile('~/tmp/git_repo/project/code/results_aurora/',name_mat));

        if( strcmp(name_model,'cppca_em') || strcmp(name_model,'cppca_affine'))
            models = ld.models_desc{idx_model}.models;
            save_mat = sprintf('%s_%s.mat',expr_desc.name,names_models{idx_m});                
            disp(save_mat)
            save(fullfile('../../results',save_mat),'models');
	    clear models;
        else        
            for idx_ETA=1:n_ETAarr        
            for idx_NV=1:n_NVarr
                NV = NVarr(idx_NV);
                Network = get_adj_graph_ex(NV);
                n_Network = numel(Network);
                for idx_Network=1:n_Network
        
                models = ld.models_desc{idx_model}.models{idx_ETA,idx_NV}{idx_Network};
                save_mat = sprintf('%s_%s_ETA%d_NV%d_%s.mat',expr_desc.name,names_models{idx_m},ETAarr(idx_ETA), NVarr(idx_NV), Network{idx_Network}.name);
                disp(save_mat)
                save(fullfile('../../results',save_mat),'models');
	        clear models;
                end
            end
            end
        end
    end
       %     ld = load( fullfile('~/tmp/git_repo/project/code/software/results_aurora/',...
        %                expr_desc.name, 
    
end
