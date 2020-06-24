close all;
clearvars;

env_setting;

%dir_summary = '~/tmp2/results20rms/summary';
%'~/tmp2/results10ssa/summary';%
%dir_summary = '~/tmp2/results10ssa/summary';
%dir_summary = '~/tmp2/results_v5uniform/summary';
%dir_summary = '~/tmp2/results_sfm_new/summary';
%dir_summary = '~/tmp2/results_primdual/summary';
%dir_summary = '~/tmp2/results_final/summary';
%dir_summary = '~/tmp2/results_final2/summary';
%dir_summary = '~/tmp2/results_supp/summary';
%dir_summary = '~/tmp2/results_supp2/summary';
dir_summary = '~/tmp2/results_supp3/summary';

names_Networks = {'Complete', 'Ring', 'Star', 'Chain', 'Cluster1', 'Cluster'};
names_Models   = {'CPPCA (svd)', 'CPPCA (em)', 'CPPCA (affine)','DPPCA','Ours (ADMM-NAP)', 'Ours (ADMM-AP)', 'F-ADMM-RS', 'ADMM-VP'};
names_expr     = {'Synthetic Gaussian', 'Synthetic Gaussian Random', 'Synthetic Cube', 'Caltech', 'Hopkins'};
    
experiments_desc = load_experiments_desc;
range_experiments = 4%[2 3 4]%2%[2 4];%[2:3];
range_models = [4 6 5]%[4 7 8 6 5];
range_Network = [1]

fontsz_title = 14;
fontsz_x     = 12;
fontsz_y     = 12;

colors = ...
[ ...
    0.9290    0.6940    0.1250; % ADMM
    0    0.4470    0.7410;% ADMM (AP)
    0.8500    0.3250    0.0980;% ADMM (NAP)
];
    

% [ ...
%     0.9290    0.6940    0.1250; % ADMM
%     0.4940    0.1840    0.5560;
%     0.4660    0.6740    0.1880;
%     0    0.4470    0.7410;% ADMM (AP)
%     0.8500    0.3250    0.0980;% ADMM (NAP)
%     %0.6350    0.0780    0.1840; ; 
%     0.3010    0.7450    0.9330;
% 
% ];

% [ ...
%     0.9290    0.6940    0.1250; % ADMM
%     0    0.4470    0.7410;% ADMM (AP)
%     0.8500    0.3250    0.0980;% ADMM (NAP)
%     %0.6350    0.0780    0.1840; ; 
%     0.4940    0.1840    0.5560;
%     0.4660    0.6740    0.1880;
%     0.3010    0.7450    0.9330;
% 
% ];

%marks = {'s', '*', '+', 'o','x'};
marks = {'s', 'o','x'};

%[ 0,0,0; ...          
%          204, 51, 0; ...
%           255, 153, 51; ...
%                          ] ./ 255;
%                      
colors_nodes = [...
 1         0         0;
 0.9290    0.6940    0.1250;
 0.4660    0.6740    0.1880;
 0    0.4470    0.7410;
 0.4940    0.1840    0.5560;
 ];

colors_nodes2 = [...
0, 51, 0
51, 204, 51
204, 204, 0
0, 153, 153
0, 102, 153 ] ./255;
                     

for idx_ex=range_experiments
    expr_desc = experiments_desc{idx_ex};
    name_expr = expr_desc.name;
    
    load(fullfile(dir_summary,sprintf('%s_summary.mat',name_expr)));
    
    NVarr  = expr_desc.NVarr;
    ETAarr = expr_desc.ETAarr;
    n_NVarr  = length(expr_desc.NVarr);
    n_ETAarr = length(expr_desc.ETAarr);
    
    if( ~exist('range_NV','var') )
        range_NV = 1:n_NVarr;
    end
    
    %n_run = expr_desc.n_run;
    n_run = 20;
    n_models = numel(expr_desc.models);

    n_params = numel(expr_desc.params);
    n_params2 = numel(expr_desc.params2);
    if n_params2==0
        n_params2 = 1;
    end
    
    for idx_ex_param=1:n_params  
    params_expr = expr_desc.params;
    params2_expr = expr_desc.params2;
        
    for idx_ex_param2=1:n_params2    
    for idx_ETA=1:n_ETAarr
        ETA = ETAarr(idx_ETA);
        
        if idx_ex > 2
            if ~exist('range_Network', 'var')
                range_Network = 1:4;
            else
                range_Network = range_Network(range_Network~=5&range_Network~=6);
            end
        else
            if ~exist('range_Network', 'var')
                range_Network = 1:6;
            end
        end
        
        for idx_Network=range_Network
                
            iter = zeros([length(range_models), n_NVarr, n_run ]);
            for i=1:length(range_models)                        
                iter(i,:,:) = squeeze(results{range_models(i),idx_ex_param,idx_ex_param2}.iter(idx_ETA, :, idx_Network, :));
            end            
            
            time = zeros([length(range_models), n_NVarr, n_run ]);
            for i=1:length(range_models)                        
                time(i,:,:) = squeeze(results{range_models(i),idx_ex_param,idx_ex_param2}.time(idx_ETA, :, idx_Network, :));
            end            
            
            err = zeros([length(range_models), n_NVarr, n_run ]);
            for i=1:length(range_models)                      
                %err(i,:,:) = squeeze(results{range_models(i),idx_ex_param,idx_ex_param2}.err_rms(idx_ETA, :, idx_Network, :));
                %err(i,:,:) = squeeze(results{range_models(i),idx_ex_param,idx_ex_param2}.err_ssa(idx_ETA, :, idx_Network, :));
                err(i,:,:) = 180/pi*squeeze(results{range_models(i),idx_ex_param,idx_ex_param2}.err_ssa(idx_ETA, :, idx_Network, :));
            end                    
                            
            idx_target = zeros( length(range_models), n_NVarr);
            err_median = zeros( length(range_models), n_NVarr);
            err_std    = zeros( length(range_models), n_NVarr);
            err_mean   = zeros( length(range_models), n_NVarr);
            ssa_median = cell(  length(range_models), n_NVarr);            
            rms_median = cell(  length(range_models), n_NVarr);                        
            obj_median = cell(  length(range_models), n_NVarr);
            eta_median = cell(  length(range_models), n_NVarr);
            iter_median= zeros( length(range_models), n_NVarr);
            iter_std   = zeros( length(range_models), n_NVarr);
            iter_mean  = zeros( length(range_models), n_NVarr);
            time_std   = zeros( length(range_models), n_NVarr);
            time_mean  = zeros( length(range_models), n_NVarr);
            for i=1:length(range_models)
                idx_model = range_models(i);
                for idx_NV=range_NV
                    
                    if( idx_ex > 2 && idx_NV > 1)
                        continue;
                    end
                    
                    if n_NVarr > 1

                    dist = abs(err(i,idx_NV,:) - median(err(i,idx_NV,:)));
                    %dist = abs(err(i,idx_NV,:) - min(err(i,idx_NV,:)));
                    %dist = abs(iter(i,idx_NV,:) - median(iter(i,idx_NV,:)));
                    [~, idx_tmp] = sort(dist,'ascend');
                    
                    idx_target(i,idx_NV)  = idx_tmp(1);

                    err_median(i,idx_NV)  = err(i,idx_NV,idx_tmp(1));
%                     err_std(i,idx_NV)     = std(err(i,idx_NV, idx_tmp(1:ceil(n_run/2))));
%                     err_mean(i,idx_NV)    = mean(err(i,idx_NV, idx_tmp(1:ceil(n_run/2))));
                    err_std(i,idx_NV)     =  std(err(i,idx_NV, :));
                    err_mean(i,idx_NV)    = mean(err(i,idx_NV, :));
                    
                    
                    iter_median(i,idx_NV) = iter(i,idx_NV,idx_tmp(1));
                    if(iter_median(i,idx_NV)>1000)
                        iter_median(i,idx_NV) = 1000;
                    end
%                     iter_std(i,idx_NV)    = std(iter(i,idx_NV, idx_tmp(1:ceil(n_run/2))));
%                     iter_mean(i,idx_NV)   = mean(iter(i,idx_NV, idx_tmp(1:ceil(n_run/2))));
                    iter_std(i,idx_NV)    =  std(iter(i,idx_NV, :));
                    iter_mean(i,idx_NV)   = mean(iter(i,idx_NV, :));
                    
                    time_std(i,idx_NV)    =  std(time(i,idx_NV, :));
                    time_mean(i,idx_NV)   = mean(time(i,idx_NV, :));

                    else
                        
                    dist = abs(err(i,:) - median(err(i,:)));
                    [~, idx_tmp] = sort(dist,'ascend');

                    idx_target(i)  = idx_tmp(1);

                    err_median(i)  = err(i,idx_tmp(1));
                    err_std(i)     =  std(err(i, idx_tmp(1:ceil(n_run*3/4))));
                    err_mean(i)    = mean(err(i, idx_tmp(1:ceil(n_run*3/4))));
%                     err_std(i)     =  std(err(i, :));
%                     err_mean(i)    = mean(err(i, :));
                                        

                    iter_median(i) = iter(i,idx_tmp(1));
                    if(iter_median(i)>1000)
                        iter_median(i) = 1000;
                    end
                    
                    iter_std(i)    =  std(iter(i, idx_tmp(1:ceil(n_run*3/4))));
                    iter_mean(i)   = mean(iter(i, idx_tmp(1:ceil(n_run*3/4))));
                    
                    time_std(i)    =  std(iter(i, idx_tmp(1:ceil(n_run*3/4))));
                    time_mean(i)   = mean(iter(i, idx_tmp(1:ceil(n_run*3/4))));

%                     iter_std(i)    =  std(iter(i, :));
%                     iter_mean(i)   = mean(iter(i, :));
%                     
%                     time_std(i)    =  std(iter(i, :));
%                     time_mean(i)   = mean(iter(i, :));

                        
                    end
                    
                    ssa_median{i,idx_NV}...
                     = results{range_models(i),idx_ex_param,idx_ex_param2}.ssa{idx_ETA, idx_NV, idx_Network, idx_tmp(1)};
                    
                    rms_median{i,idx_NV}...
                     = results{range_models(i),idx_ex_param,idx_ex_param2}.rms{idx_ETA, idx_NV, idx_Network, idx_tmp(1)};
                 
                    obj_median{i,idx_NV}...
                     = results{range_models(i),idx_ex_param,idx_ex_param2}.obj{idx_ETA, idx_NV, idx_Network, idx_tmp(1)};                 
                 
                    eta_median{i,idx_NV}...
                     = results{range_models(i),idx_ex_param,idx_ex_param2}.eta{idx_ETA, idx_NV, idx_Network, idx_tmp(1)};                 
                end
            end
if 1
            name_target = 'Mean of Subspace Error';
            name_xlabel = 'number of Nodes';
            name_ylabel = 'subspace error';
   
            errorbar_groups(err_mean, err_std, ...
                            'bar_colors', colors(1:length(range_models),:), ...
                            'errorbar_colors', colors(1:length(range_models),:)*0.5,...
                            'bar_names', cellstr(num2str(NVarr')));
                        
            if( idx_ex < 4 )
            str_title = sprintf('%s %s \n(ETA:%d %s)',names_expr{idx_ex},name_target,ETA,names_Networks{idx_Network});
            else
            str_title = sprintf('%s %s \n(ETA:%d Obj:%d %s)',names_expr{idx_ex},name_target,ETA,params_expr{idx_ex_param}{1},names_Networks{idx_Network});                    
            end
            title( str_title, 'FontSize', fontsz_title);
            xlabel(name_xlabel,'FontSize',fontsz_x);
            ylabel(name_ylabel,'FontSize',fontsz_y);
            %legend( {names_Models{range_models}}, 'Location','NE');
            %ylim([0 max(max(err_mean + err_std))*1.1]);
            %ylim([0 max(0.5, max(max(err_mean + err_std))*1.1)]);
            ylim([0 1]);
            
            if( idx_ex < 4 )
            str_title = sprintf('%s %s (ETA:%d %s)',names_expr{idx_ex},name_target,ETA,names_Networks{idx_Network});
            else
            str_title = sprintf('%s %s (ETA:%d Obj:%d %s)',names_expr{idx_ex},name_target,ETA,params_expr{idx_ex_param}{1},names_Networks{idx_Network});                    
            end
            file_img = strrep(str_title, ' ', '_');
            file_img = strrep(file_img, ':', '');
            file_img = strrep(file_img, '(', '');                    
            file_img = strrep(file_img, ')', '');                                        
            saveas(gcf, sprintf('../../img/%s',file_img), 'epsc');
            saveas(gcf, sprintf('../../img/%s',file_img), 'png');

            close gcf;
            
            name_target = 'Mean of Number of Iterations';
            name_xlabel = 'number of Nodes';
            name_ylabel = 'number of iterations';

            errorbar_groups(iter_mean, iter_std, ...
                'bar_colors', colors(1:length(range_models),:), ...
                'errorbar_colors', colors(1:length(range_models),:)*0.5,...
                'bar_names', cellstr(num2str(NVarr')));

            if( idx_ex < 4 )
            str_title = sprintf('%s %s \n(ETA:%d %s)',names_expr{idx_ex},name_target,ETA,names_Networks{idx_Network});
            else
            str_title = sprintf('%s %s \n(ETA:%d Obj:%d %s)',names_expr{idx_ex},name_target,ETA,params_expr{idx_ex_param}{1},names_Networks{idx_Network});                    
            end
            title( str_title, 'FontSize', fontsz_title);
            xlabel(name_xlabel,'FontSize',fontsz_x);
            ylabel(name_ylabel,'FontSize',fontsz_y);
            %legend( {names_Models{range_models}}, 'Location','SW');
            ylim([0 max(max(iter_mean + iter_std))*1.1]);
            %ylim([0 max(120, max(max(iter_mean + iter_std))*1.1)]);
            
            
            if( idx_ex < 4 )
            str_title = sprintf('%s %s (ETA:%d %s)',names_expr{idx_ex},name_target,ETA,names_Networks{idx_Network});
            else
            str_title = sprintf('%s %s (ETA:%d Obj:%d %s)',names_expr{idx_ex},name_target,ETA,params_expr{idx_ex_param}{1},names_Networks{idx_Network});                    
            end
            file_img = strrep(str_title, ' ', '_');
            file_img = strrep(file_img, ':', '');
            file_img = strrep(file_img, '(', '');                    
            file_img = strrep(file_img, ')', '');                                        
            saveas(gcf, sprintf('../../img/%s',file_img), 'epsc');
            saveas(gcf, sprintf('../../img/%s',file_img), 'png');

            close gcf;
            
            name_target = 'Mean of Running Time';
            name_xlabel = 'number of Nodes';
            name_ylabel = 'Running Time (sec)';

            errorbar_groups(time_mean, time_std, ...
                'bar_colors', colors(1:length(range_models),:), ...
                'errorbar_colors', colors(1:length(range_models),:)*0.5,...
                'bar_names', cellstr(num2str(NVarr')));

            if( idx_ex < 4 )
            str_title = sprintf('%s %s \n(ETA:%d %s)',names_expr{idx_ex},name_target,ETA,names_Networks{idx_Network});
            else
            str_title = sprintf('%s %s \n(ETA:%d Obj:%d %s)',names_expr{idx_ex},name_target,ETA,params_expr{idx_ex_param}{1},names_Networks{idx_Network});                    
            end
            title( str_title, 'FontSize', fontsz_title);
            xlabel(name_xlabel,'FontSize',fontsz_x);
            ylabel(name_ylabel,'FontSize',fontsz_y);
            %legend( {names_Models{range_models}}, 'Location','SW');
            ylim([0 max(max(iter_mean + iter_std))*1.1]);
            %ylim([0 max(120, max(max(iter_mean + iter_std))*1.1)]);
            %ylim([0 inf]);
            
            if( idx_ex < 4 )
            str_title = sprintf('%s %s (ETA:%d %s)',names_expr{idx_ex},name_target,ETA,names_Networks{idx_Network});
            else
            str_title = sprintf('%s %s (ETA:%d Obj:%d %s)',names_expr{idx_ex},name_target,ETA,params_expr{idx_ex_param}{1},names_Networks{idx_Network});                    
            end
            file_img = strrep(str_title, ' ', '_');
            file_img = strrep(file_img, ':', '');
            file_img = strrep(file_img, '(', '');                    
            file_img = strrep(file_img, ')', '');                                        
            saveas(gcf, sprintf('../../img/%s',file_img), 'epsc');
            saveas(gcf, sprintf('../../img/%s',file_img), 'png');

            close gcf;

end
if 0
            for idx_NV=range_NV
                
                if( idx_ex > 2 && idx_NV > 1)
                    continue;
                end
                
                NV=NVarr(idx_NV);
                Networks = get_adj_graph_ex(NV);
                n_Networks = length(Networks);        

                h1 = figure();
                
                haxes = axes;
                hold on;            
                for i=1:length(range_models)
%                     plot( 1:iter_median(i,idx_NV), ...
%                           rms_median{i,idx_NV}(1:iter_median(i,idx_NV)), ...
%                           [marks{i} '-'], 'color', colors(i,:), 'LineWidth', 2 );
                    
                    plot( 1:iter_median(i,idx_NV), ...
                          180/pi*ssa_median{i,idx_NV}(1:iter_median(i,idx_NV)), ...
                          [marks{i} '-'], 'color', colors(i,:), 'LineWidth', 2 );
                end
                
                %set(haxes,'XScale','log');
                set(haxes,'YScale','log');
                
                xlim([0, max(iter_median(:,idx_NV))]);
                ylim([0, inf]);
                
                hold off;
                
                name_target = 'Subspace Error';
                name_xlabel = 'iteration';
                name_ylabel = 'subspace error';

                if( idx_ex < 4 )
                str_title = sprintf('%s %s \n(ETA:%d Nodes:%d %s)',names_expr{idx_ex},name_target,ETA,NV,names_Networks{idx_Network});
                else
                str_title = sprintf('%s %s \n(ETA:%d Obj:%d %s)',names_expr{idx_ex},name_target,ETA,params_expr{idx_ex_param}{1},names_Networks{idx_Network});                    
                end
                title( str_title, 'FontSize', fontsz_title);
                xlabel(name_xlabel,'FontSize',fontsz_x);
                ylabel(name_ylabel,'FontSize',fontsz_y);
                legend( {names_Models{range_models}}, 'Location','SW');
                
                if( idx_ex < 4 )
                str_title = sprintf('%s %s (ETA:%d Nodes:%d %s)',names_expr{idx_ex},name_target,ETA,NV,names_Networks{idx_Network});
                else
                str_title = sprintf('%s %s (ETA:%d Obj:%d %s)',names_expr{idx_ex},name_target,ETA,params_expr{idx_ex_param}{1},names_Networks{idx_Network});                    
                end
                file_img = strrep(str_title, ' ', '_');
                file_img = strrep(file_img, ':', '');
                file_img = strrep(file_img, '(', '');                    
                file_img = strrep(file_img, ')', '');                                        
                saveas(gcf, sprintf('../../img/%s',file_img), 'epsc');
                saveas(gcf, sprintf('../../img/%s',file_img), 'png');
                
                close gcf;
                
                
                
%                 fprintf('Model: AP\n');                
%                 eta_tmp = eta_median{end-1,idx_NV}(1:iter_median(end-1,idx_NV),:,:);
%                 eta_sum = sum(sum(sum(eta_tmp)));
%                 n_edges = sum(sum(sum(0<eta_tmp)));
%                 fprintf('ETA: 5~20\n');
%                 fprintf('NV\t%d\tNet\t%s\titer\t%d\tn_Edges\t%d\tETA\t%d\twEta\t%f\n',NV,Networks{idx_Network}.name, iter_median(end,idx_NV), n_edges, ETA, eta_sum/n_edges);
% 
%                 eta_tmp(eta_tmp>ETA) = (eta_tmp(eta_tmp>ETA) - ETA)/10 * 5 + ETA;
%                 eta_sum = sum(sum(sum(eta_tmp)));
%                 n_edges = sum(sum(sum(0<eta_tmp)));
%                 fprintf('ETA: 5~15\n');                
%                 fprintf('NV\t%d\tNet\t%s\titer\t%d\tn_Edges\t%d\tETA\t%d\twEta\t%f\n',NV,Networks{idx_Network}.name, iter_median(end,idx_NV), n_edges, ETA, eta_sum/n_edges);
                 
%                 fprintf('Model: NAP\n');                
%                 eta_tmp = eta_median{end,idx_NV}(1:iter_median(end,idx_NV),:,:);
%                 eta_sum = sum(sum(sum(eta_tmp)));
%                 n_edges = sum(sum(sum(0<eta_tmp)));
%                 fprintf('ETA: 5~20\n');
%                 fprintf('NV\t%d\tNet\t%s\titer\t%d\tn_Edges\t%d\tETA\t%d\twEta\t%f\n',NV,Networks{idx_Network}.name, iter_median(end,idx_NV), n_edges, ETA, eta_sum/n_edges);
% 
%                 eta_tmp(eta_tmp>ETA) = (eta_tmp(eta_tmp>ETA) - ETA)/10 * 5 + ETA;
%                 eta_sum = sum(sum(sum(eta_tmp)));
%                 n_edges = sum(sum(sum(0<eta_tmp)));
%                 fprintf('ETA: 5~15\n');                
%                 fprintf('NV\t%d\tNet\t%s\titer\t%d\tn_Edges\t%d\tETA\t%d\twEta\t%f\n',NV,Networks{idx_Network}.name, iter_median(end,idx_NV), n_edges, ETA, eta_sum/n_edges);
                
                
                if( NV==6 )
                    for id = [2 3]
                        
                    if( id == 2 )
                    name_target = ['Changes of ETA' ' AP'] ;
                    else
                    name_target = ['Changes of ETA' ' NAP'] ;                        
                    end
                    name_xlabel = 'iteration';
                    name_ylabel = 'penalty constraint (eta)';
                    %iters = [5 25 40 50 55];
                     %iters = [5 15 35 40 44]%5:10:iter_median(id,idx_NV);
                     %iters = [ 5 25 40 50 55];
                     iters = [5 20 45 46 55];
                     
                     iters = iters(iters < iter_median(id,idx_NV));                     
                     n_iters = length(iters);
%                     hold on;
%                     ah = axes('Position',[.15 .6 .2 .2]); 
%                     quiver(ah, 0.75, 0.75, 1, 1, 'k', 'filled', 'MaxHeadSize',2, 'LineWidth',2);
%                     set(ah, 'Box', 'off');
%                     set(ah, 'XTick',[]);
%                     set(ah, 'YTick',[]);
%                     set(ah,'visible','off');
%                     axis([0 2 0 2]);
%                     hold on;
%                     show_topology( NV, 5, ETA, eta_median{end,idx_NV}, colors_nodes );

                

                    figure();
                    if( idx_ex < 4 )
                        str_title = sprintf('%s %s \n(ETA:%d Nodes:%d %s)',names_expr{idx_ex},name_target,ETA,NV,names_Networks{idx_Network});
                    else
                        str_title = sprintf('%s %s \n(ETA:%d Obj:%d %s)',names_expr{idx_ex},name_target,ETA,params_expr{idx_ex_param}{1},names_Networks{idx_Network});                    
                    end     
                    iter_tmp = 1:iter_median(id,idx_NV);
                    for a=1:NV
                        subplot(NV+ceil(n_iters/6),1,a+ceil(n_iters/6));
                        hold on;
                        for b=1:NV
                            if a==b 
                                continue;
                            end                  
                            if eta_median{id,idx_NV}(iter_tmp,a,b)==0 
                                continue;
                            end
                            stairs( iter_tmp - 0.5, ...
                                  eta_median{id,idx_NV}(iter_tmp,a,b),...
                                  'color', colors_nodes(b,:), 'LineWidth', 1.5 );   
                        end

                        xlim([0.5, iter_median(id,idx_NV)-0.5]);
                        %ylim([ETA*0.5-1 ETA*2+1]);
                        ylim([ETA*0.5-1 max(max(max(eta_median{id,idx_NV}))) + 1]);
                        if a==1
                            %title( str_title, 'FontSize', fontsz_title);                    
                        end
                        if a==ceil(NV/2)
                            ylabel(name_ylabel,'FontSize',fontsz_y);
                        end
                        if a==NV
                            xlabel(name_xlabel,'FontSize',fontsz_x);
                        end
                        if a < NV
                            set(gca,'XTick',[]);
                        else

                        end                    
                    end
                 
                    n_iters = length(iters);
                    for idx_iter=1:n_iters
                        iter = iters(idx_iter);
                        if( n_iters <= 6 )
                            subplot(n_iters,n_iters,idx_iter);
                        elseif( n_iters <=12 )
                            subplot(ceil(n_iters/2)+3,ceil(n_iters/2),idx_iter);
                        else
                            subplot(ceil(n_iters/3)+5,ceil(n_iters/3),idx_iter);
                        end

                        hold on;

                        show_topology( NV, iter, ETA, eta_median{id,idx_NV}, colors_nodes );
                    end
                    
                    ah=axes('visible','off');                                        
                    hold on;
                    for p=1:n_iters                        
                        line('XData',[iters(p) iters(p)],'YData',[0 (1-ceil(n_iters/6)*0.12)-0.05], 'parent',ah,'linewidth',1,'LineStyle','-.','color',[1 0 0]);%[0 0.4470 0.7410]
                    end
                    
                    for a=1:NV
                    h = plot(ah,0.75,(1-ceil(n_iters/6)*0.12)*(1-(a)/(NV))+ 0.06,'o','LineWidth',4,'MarkerEdgeColor',colors_nodes(a,:),'MarkerFaceColor',colors_nodes(a,:));
                        text(   1.2, (1-ceil(n_iters/6)*0.12)*(1-(a)/(NV))+0.075,num2str(a),'color',colors_nodes(a,:), 'FontSize',12,'FontWeight','bold');
                    end
                    xlim(ah,[0.5, iter_median(id,idx_NV)-0.5]);
                    ylim(ah,[0 1]);
                    
                    hold off;
                    if( idx_ex < 4 )
                    str_title = sprintf('%s %s (ETA:%d Nodes:%d %s)',names_expr{idx_ex},name_target,ETA,NV,names_Networks{idx_Network});
                    else
                    str_title = sprintf('%s %s (ETA:%d Obj:%d %s)',names_expr{idx_ex},name_target,ETA,params_expr{idx_ex_param}{1},names_Networks{idx_Network});                    
                    end                
                    file_img = strrep(str_title, ' ', '_');
                    file_img = strrep(file_img, ':', '');
                    file_img = strrep(file_img, '(', '');                    
                    file_img = strrep(file_img, ')', '');                                        
                    saveas(gcf, sprintf('../../img/%s',file_img), 'epsc');
                    saveas(gcf, sprintf('../../img/%s',file_img), 'png');
                
                    close gcf;
                    end
                end                
                
            end
end            
        end
    end
    end
    end
end

% idx_target = zeros( length(range_models), n_NVarr);
%                 err_median = zeros( length(range_models), n_NVarr);
%                 iter_median= zeros( length(range_models), n_NVarr);
%                 for idx_model=1:length(range_models)
%                     for idx_NV=1:n_NVarr
%                         [~, idx_tmp] ...
%                          = min(abs(err(idx_model,idx_NV,:) - median(err(idx_model,idx_NV,:))));
%                      
%                         idx_target(idx_model,idx_NV)  = idx_tmp;
%                         err_median(idx_model,idx_NV)  = err(idx_model,idx_NV,idx_tmp);
%                         iter_median(idx_model,idx_NV) = iter(idx_model,idx_NV,idx_tmp);
%                     end
%                 end
%  errorbar_groups(iter_median, iter_std, ...
%                 'bar_colors', colors(1:length(range_models),:), ...
%                 'errorbar_colors', colors(1:length(range_models),:)*0.5 );
% name_target = 'Number of Iterations';
% name_xlabel = '# of nodes';
% name_ylabel = '# of iterations';
% ylim([0 inf]);
% set(gca,'XTickLabel', cellstr(num2str(NVarr(:))) );
% str_title = sprintf('%s %s \n(ETA:%d, %s)',names_expr{idx_ex},name_target,ETA,names_Networks{idx_Network});
% title( str_title, 'FontSize', fontsz_title);
% xlabel(name_xlabel,'FontSize',fontsz_x);
% ylabel(name_ylabel,'FontSize',fontsz_y);
% legend( {names_Models{range_models}}, 'Location','NW');
% saveas(gcf, sprintf('../../img/%s',str_title), 'png');
% hold off;
% 
% err_std     = std(err,[],3);
% errorbar_groups(err_median, err_std, ...
%                 'bar_colors', colors(1:length(range_models),:), ...
%                 'errorbar_colors', colors(1:length(range_models),:)*0.5 );
% name_target = 'Subspace Error';
% name_xlabel = '# of nodes';
% name_ylabel = 'subspace error';
% ylim([0 inf]);                                
% set(gca,'XTickLabel', cellstr(num2str(NVarr(:))) );
% str_title = sprintf('%s %s \n(ETA:%d, %s)',names_expr{idx_ex},name_target,ETA,names_Networks{idx_Network});
% title( str_title, 'FontSize', fontsz_title);
% xlabel(name_xlabel,'FontSize',fontsz_x);
% ylabel(name_ylabel,'FontSize',fontsz_y);
% legend( {names_Models{range_models}}, 'Location','NW');
% saveas(gcf, sprintf('../../img/%s',str_title), 'png');
% hold off;


%                 figure();
%                 hold on;
%                 x_step_size = ceil(length(range_models)/3);
%                 for i=1:length(range_models)
%                     idx_model = range_models(i);
% 
%                     position = x_step_size*(1:n_NVarr) + ...
%                                (0.3*(-(length(range_models)-1)/2)) + ...
%                                0.3*(i-1);
%                     boxplot(squeeze(err(i,:,:))',...                        
%                         'labels',  cellstr(num2str(NVarr(:))), ...
%                         'colors',colors(i,:)*0.8,'symbol','g+','Position',position,'widths',0.2); 
%                     %set(gca,'XTickLabel',{' '});
%                 end
%                 name_target = 'Subspace Error';
%                 name_xlabel = '# of nodes';
%                 name_ylabel = 'subspace error';
%                 ylim([0 inf]);                                
%                 set(gca,'XTickLabel', cellstr(num2str(NVarr(:))) );
%                 str_title = sprintf('%s %s \n(ETA:%d, %s)',names_expr{idx_ex},name_target,ETA,names_Networks{idx_Network});
%                 title( str_title, 'FontSize', fontsz_title);
%                 xlabel(name_xlabel,'FontSize',fontsz_x);
%                 ylabel(name_ylabel,'FontSize',fontsz_y);
%                 legend( {names_Models{range_models}}, 'Location','NW');
%                 saveas(gcf, sprintf('../../img/%s',str_title), 'png');
%                 hold off;
