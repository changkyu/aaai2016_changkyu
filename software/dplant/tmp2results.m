load('./supp_caltech.mat','-ASCII');

dir_result = '~/tmp2/results_supp4';

NV = 5;
ETA = 10;
n_ETAarr = 1;
n_NVarr = 1;
n_run = 20;

results = cell(8,5,1);

name_expr = 'ex08_sfm_caltech';

for idx_model=4:6
    
    if( idx_model==4 )
        %filepath_mat = '../../../results_final/dppca.mat';
        %ld = load(filepath_mat);
        %cm = ld.cm1;
        cm = cm1;
        name_model = 'dppca';
    elseif( idx_model==5 )
        %filepath_mat = '../../../results_final/dppca_ant.mat';
        %ld = load(filepath_mat);
        %cm = ld.cm1;
        cm = cm3;
        name_model = 'dppca_ant';

    elseif( idx_model==6 )
%         filepath_mat = '../../../results_final/dppca_ant_Tmax.mat';
%         ld = load(filepath_mat);
%         cm = ld.cm1;
        cm = cm2;        
        name_model = 'dppca_ant_Tmax';
    end
    
    
    for idx_ex_param=1:5
        idx_ex_param2 = 1;
        
        results{idx_model,idx_ex_param,idx_ex_param2}.iter = zeros(n_ETAarr,n_NVarr,6,n_run);
        results{idx_model,idx_ex_param,idx_ex_param2}.time = zeros(n_ETAarr,n_NVarr,6,n_run);
        %results{idx_model,idx_ex_param,idx_ex_param2}.err  = zeros(n_ETAarr,n_NVarr,6,n_run);
        results{idx_model,idx_ex_param,idx_ex_param2}.err_rms  = zeros(n_ETAarr,n_NVarr,6,n_run);
        results{idx_model,idx_ex_param,idx_ex_param2}.err_ssa  = zeros(n_ETAarr,n_NVarr,6,n_run);
        results{idx_model,idx_ex_param,idx_ex_param2}.obj  =  cell(n_ETAarr,n_NVarr,6,n_run);            
        results{idx_model,idx_ex_param,idx_ex_param2}.ssa  =  cell(n_ETAarr,n_NVarr,6,n_run);            
        results{idx_model,idx_ex_param,idx_ex_param2}.rms  =  cell(n_ETAarr,n_NVarr,6,n_run);            
        results{idx_model,idx_ex_param,idx_ex_param2}.eta  =  cell(n_ETAarr,n_NVarr,6,n_run);            
        
        idx_ETA=1;
        idx_NV =1;
        for idx_Network=1:4
            for idx_run=1:n_run
                model = cm{idx_ex_param,idx_Network,idx_run};
                
                err_ssa = model.ssaArray(end);
                err_rms = model.rmsArray(end);

                results{idx_model,idx_ex_param,idx_ex_param2}.iter(idx_ETA,idx_NV,idx_Network,idx_run) = model.eITER;
                results{idx_model,idx_ex_param,idx_ex_param2}.time(idx_ETA,idx_NV,idx_Network,idx_run) = model.eTIME;
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
            end
        end                
    end
end

if( ~exist( fullfile(dir_result, 'summary'), 'dir' ) )
    mkdir(fullfile(dir_result, 'summary'));
end
save(fullfile(dir_result, 'summary', sprintf('%s_summary.mat',name_expr)), 'results','-v7.3');
clear results;    
