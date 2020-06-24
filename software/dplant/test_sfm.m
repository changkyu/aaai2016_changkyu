clear; close all;

% Choose random seed: optional setting to reproduce numbers.
s = RandStream('mt19937ar','Seed',1);
RandStream.setGlobalStream(s);
reset(s,0);

%dataset = 'caltech';
dataset = 'hopkins';

range_Network = 1:4;
if( strcmp(dataset,'caltech') )    
%    range_obj = 1:5;
    range_obj = 1:5;
    n_run = 1;
elseif(strcmp(dataset, 'hopkins'))
    range_obj = 1:135;
    n_run = 1;
end

vis = false;
conv_counter = 4;
conv_counter_rms = 4;

if( vis )    
    range_Network = 2;
    range_obj = 3;
    n_run = 1;
    filepath_save_mat = sprintf('./draw_%s_obj%d.mat', dataset, range_obj);
else
    filepath_save_mat = ['./final_' dataset '.mat']; 
end

cm0 = cell(5,6,n_run);
cm1 = cell(5,6,n_run);
cm2 = cell(5,6,n_run);
cm3 = cell(5,6,n_run);
cm4 = cell(5,6,n_run);
cm5 = cell(5,6,n_run);
cm6 = cell(5,6,n_run);
% cm7 = cell(5,6,n_run);
% cm8 = cell(5,6,n_run);
GT  = cell(5,6,n_run);

for jj = range_Network
    for ii = range_obj
        load(['D:/Work/NIPS2015/project/code/dataset/' dataset '/' sprintf('%03d',ii) '.mat']);
        X = measurement_matrix;
        mm_mat_t = X';
        [~, N] = size(X);

        centroid = mean(X, 2);
        X_trans = X - repmat(centroid, [1, N]);
        X = X_trans';

        [D, N] = size(X);

        M = 3;

        % ETA: Learning rate
        ETA = 10;

        % Node assignment to each sample
        NV = 5;
        V = get_sample_assign(NV, N);

        for r=1:n_run
            % E: Adjacency graph topology (1:complete, 2:ring, 3:star, 4:chain)
            GRAPHS = get_adj_graph(NV);
            E = GRAPHS{jj};

            % For structure from motion applications only, better to enforce zero mean
            ZeroMean = true;

            % Convergence thresholds (typically you don't have to change these)
            THRESHc = 1e-3;     % threshold (objective) for PPCA
            THRESHd = 1e-3;     % threshold (objective) for D-PPCA
            THRESHda = 1e-4;    % threshold (absolute)  for D-PPCA
            THRESHdr = 1e-2;    % threshold (relative)  for D-PPCA

            % Show progress per this number of iterations
            objfreq_c = 0;      % show objective per this iterations (PPCA)
            objfreq_d = 10;    % show objective per this iterations (D-PPCA)

            [sfm_M, sfm_S, sfm_b, U3, W3, V3] = sfm_affine(X');
            GT{ii,jj,r} = sfm_S';

            % centralized model initialization: 
            %   a random projection, zero mean and a small variance
            m_init_c = get_init_value(X, M, 'ModelType', 'cppca');

            % distributed model initialization:
            %   initializes with a local PPCA result
            m_init_d = get_init_value(X, M, 'ModelType', 'sfm_d', ...
                'NumberOfNodes', NV, 'SampleAssignVec', V);

            % Run PPCA
            disp('* PPCA (EM) *');
            cm0{ii,jj,r} = cppca_em(X, M, 'InitModel', m_init_c, 'ShowObjPer', objfreq_c);

            % Run D-PPCA
            disp('* D-PPCA *');
            cm1{ii,jj,r} = dppca(X, M, V, E, 'Eta', ETA, ...
                'InitModel', m_init_d, 'Threshold', THRESHd, 'ShowObjPer', objfreq_d, ...
                'ThresholdA', THRESHda, 'ThresholdR', THRESHdr, ...
                'ZeroMean', ZeroMean, 'MaxIter', 100, ...
                'UseResidual', false, 'W_GT', sfm_S', ...
                'ConvCounter', conv_counter, 'ConvCounterRMS', conv_counter_rms);

            % Run D-PPCA
            disp('* D-PPCA AP  *');
            cm2{ii,jj,r} = dppca_ant(X, M, V, E, 'Eta', ETA, ...
                'InitModel', m_init_d, 'Threshold', THRESHd, 'ShowObjPer', objfreq_d, ...    
                'UseResidual', false,...
                'ZeroMean', ZeroMean, 'MaxIter', 100, 'W_GT', sfm_S','T',50, ...
                'ModelType', 'dppca_AP_homo', ...
                'ConvCounter', conv_counter, 'ConvCounterRMS', conv_counter_rms);

            % Run D-PPCA
            disp('* D-PPCA NAP *');
            cm3{ii,jj,r} = dppca_ant(X, M, V, E, 'Eta', ETA, ...
                'InitModel', m_init_d, 'Threshold', THRESHd, 'ShowObjPer', objfreq_d, ...    
                'UseResidual', false,...
                'ZeroMean', ZeroMean, 'MaxIter', 100, 'W_GT', sfm_S', 'vis', vis, ...
                'ModelType', 'dppca_NAP_homo', ...
                'ConvCounter', conv_counter, 'ConvCounterRMS', conv_counter_rms);

%             % Run D-PPCA
%             disp('* D-PPCA VP *');
%             cm4{ii,jj,r} = dppca_ant(X, M, V, E, 'Eta', ETA, ...
%                 'InitModel', m_init_d, 'Threshold', THRESHd, 'ShowObjPer', objfreq_d, ...    
%                 'UseResidual', true, 'T', 50, ...
%                 'ZeroMean', ZeroMean, 'MaxIter', 100, 'W_GT', sfm_S', 'vis', vis, ...
%                 'ModelType', 'dppca_localVP_hetero', 'VaryETA_mu', 0.1, 'VaryETA_tau', 1, ...
%                 'ThresholdA', THRESHda, 'ThresholdR', THRESHdr, 'ConvCounter', conv_counter, 'ConvCounterRMS', conv_counter_rms);

            % Run D-PPCA
            disp('* D-PPCA VP *');
            cm4{ii,jj,r} = dppca_ant(X, M, V, E, 'Eta', ETA, ...
                'InitModel', m_init_d, 'Threshold', THRESHd, 'ShowObjPer', objfreq_d, ...    
                'UseResidual', false, 'T', 50, ...
                'ZeroMean', ZeroMean, 'MaxIter', 100, 'W_GT', sfm_S', 'vis', vis, ...
                'ModelType', 'dppca_localVP_homo', 'VaryETA_mu', 0.1, 'VaryETA_tau', 1, ...
                'ThresholdA', THRESHda, 'ThresholdR', THRESHdr, ...
                'ConvCounter', conv_counter, 'ConvCounterRMS', conv_counter_rms);

            % Run D-PPCA
            disp('* D-PPCA VP-AP *');
            cm5{ii,jj,r} = dppca_ant(X, M, V, E, 'Eta', ETA, ...
                'InitModel', m_init_d, 'Threshold', THRESHd, 'ShowObjPer', objfreq_d, ...    
                'UseResidual', false, 'T', 50, ...
                'ZeroMean', ZeroMean, 'MaxIter', 100, 'W_GT', sfm_S', 'vis', vis, ...
                'ModelType', 'dppca_localVPAP_homo', 'VaryETA_mu', 0.1, 'VaryETA_tau', 1, ...
                'ThresholdA', THRESHda, 'ThresholdR', THRESHdr, ...
                'ConvCounter', conv_counter, 'ConvCounterRMS', conv_counter_rms);

            % Run D-PPCA
            disp('* D-PPCA VP-NAP *');
            cm6{ii,jj,r} = dppca_ant(X, M, V, E, 'Eta', ETA, ...
                'InitModel', m_init_d, 'Threshold', THRESHd, 'ShowObjPer', objfreq_d, ...    
                'UseResidual', false, ...
                'ZeroMean', ZeroMean, 'MaxIter', 100, 'W_GT', sfm_S', 'vis', vis, ...
                'ModelType', 'dppca_localVPNAP_homo', 'VaryETA_mu', 0.1, 'VaryETA_tau', 1, ...
                'ThresholdA', THRESHda, 'ThresholdR', THRESHdr, ...
                'ConvCounter', conv_counter, 'ConvCounterRMS', conv_counter_rms);
            
%             % Run D-PPCA
%             disp('* D-PPCA FRS (global) *');
%             cm7{ii,jj,r} = dppca_ant(X, M, V, E, 'Eta', ETA, ...
%                'InitModel', m_init_d, 'Threshold', THRESHd, 'ShowObjPer', objfreq_d, ...    
%                'UseResidual', false, ...
%                'ZeroMean', ZeroMean, 'MaxIter', 100, 'W_GT', sfm_S', ...
%                'ModelType', 'dppca_globalFRS', 'ConvCounter', conv_counter, 'ConvCounterRMS', conv_counter_rms);            
%             
%             % Run D-PPCA
%             disp('* D-PPCA FRS (local) *');
%             cm8{ii,jj,r} = dppca_ant(X, M, V, E, 'Eta', ETA, ...
%                'InitModel', m_init_d, 'Threshold', THRESHd, 'ShowObjPer', objfreq_d, ...    
%                'UseResidual', false, ...
%                'ZeroMean', ZeroMean, 'MaxIter', 100, 'W_GT', sfm_S', ...
%                'ModelType', 'dppca_localFRS', 'ConvCounter', conv_counter, 'ConvCounterRMS', conv_counter_rms);            

            fprintf('%d %d %d\n',ii,jj,r);
        end
    end
end

save(filepath_save_mat,'GT','cm0','cm1','cm2','cm3', 'cm4','cm5','cm6','-v7.3');

graph_names = {'complete', 'ring', 'star', 'chain'};
object_names = cell(max(range_obj),1);
if( strcmp(dataset,'caltech') )    
    object_names ={'BallSander', 'BoxStuff', 'Rooster', 'Standing', 'StorageBin'};
elseif(strcmp(dataset, 'hopkins'))
    for idx = 1:length(range_obj)
        object_names{idx} = sprintf('%03d', idx);
    end
end
method_names = {'ADMM', 'ADMM-AP', 'ADMM-NAP', ...
    'ADMM-VP', 'ADMM-VP+AP', 'ADMM-VP+NAP'};%, ...
    %'Fast ADMM w/ Restart (global)', 'Fast ADMM w/ Restart (local)'};
% subspace
for r = 1
    figure;
    for jj = range_Network
        for idx = 1:length(range_obj)
            ii=range_obj(idx);
            %subplot(1,length(range_obj),(jj-1)*length(range_obj)+idx);
            clf;
            hold on;
            loglog(rad2deg(cm1{ii,jj,r}.ssaArray), 'g');
            loglog(rad2deg(cm2{ii,jj,r}.ssaArray), 'b');
            loglog(rad2deg(cm3{ii,jj,r}.ssaArray), 'r');    
            loglog(rad2deg(cm4{ii,jj,r}.ssaArray), 'c');
            loglog(rad2deg(cm5{ii,jj,r}.ssaArray), 'b--');
            loglog(rad2deg(cm6{ii,jj,r}.ssaArray), 'r--');    
%             loglog(rad2deg(cm7{ii,jj,r}.ssaArray), 'k-');
%             loglog(rad2deg(cm8{ii,jj,r}.ssaArray), 'k--');
            set(gca,'xscale','log');
            set(gca,'yscale','log');
            ylim([0, 100]);
            grid on;

            title(sprintf('%s / %s', object_names{ii}, graph_names{jj}));
            xlabel('iterations');
            ylabel('maximum subspace angle (degrees)');
            legend(method_names, 'Location', 'SouthWest');
            
            saveas(gcf, sprintf('fig/eps/%s_Subspace_%s_%s_%02d.eps', ...
                dataset, object_names{idx}, graph_names{jj}, r), 'epsc');
            saveas(gcf, sprintf('fig/png/%s_Subspace_%s_%s_%02d.png', ...
                dataset, object_names{idx}, graph_names{jj}, r), 'png');
        end
    end
    if ~vis, close gcf; end
end

% reprojection error
for r = 1
    figure;
    for jj = range_Network
        for idx = 1:length(range_obj)
            ii=range_obj(idx);
            %subplot(1,length(range_obj),(jj-1)*length(range_obj)+idx);
            clf;
            hold on;
            loglog((cm1{ii,jj,r}.rmsArray)./(N*D), 'g');
            loglog((cm2{ii,jj,r}.rmsArray)./(N*D), 'b');
            loglog((cm3{ii,jj,r}.rmsArray)./(N*D), 'r');    
            loglog((cm4{ii,jj,r}.rmsArray)./(N*D), 'c');
            loglog((cm5{ii,jj,r}.rmsArray)./(N*D), 'b--');
            loglog((cm6{ii,jj,r}.rmsArray)./(N*D), 'r--');    
%             loglog(rad2deg(cm7{ii,jj,r}.rmsArray), 'k-');
%             loglog(rad2deg(cm8{ii,jj,r}.rmsArray), 'k--');
            set(gca,'xscale','log');
            set(gca,'yscale','log');
            grid on;

            title(sprintf('%s / %s', object_names{ii}, graph_names{jj}));
            xlabel('iterations');
            ylabel('reprojection error');
            legend(method_names, 'Location', 'SouthWest');
            
            saveas(gcf, sprintf('fig/eps/%s_ReprojectionError_%s_%s_%02d.eps', ...
                dataset, object_names{idx}, graph_names{jj}, r), 'epsc');
            saveas(gcf, sprintf('fig/png/%s_ReprojectionError_%s_%s_%02d.png', ...
                dataset, object_names{idx}, graph_names{jj}, r), 'png');
        end
    end
    if ~vis, close gcf; end
end

% objective
for r = 1
    figure;
    for jj = range_Network
        for idx = 1:length(range_obj)
            ii=range_obj(idx);
            %subplot(1,length(range_obj),(jj-1)*length(range_obj)+idx);
            clf;
            hold on;
            semilogx((cm1{ii,jj,r}.objArray(:,end)), 'g');
            semilogx((cm2{ii,jj,r}.objArray(:,end)), 'b');
            semilogx((cm3{ii,jj,r}.objArray(:,end)), 'r');    
            semilogx((cm4{ii,jj,r}.objArray(:,end)), 'c');
            semilogx((cm5{ii,jj,r}.objArray(:,end)), 'b--');
            semilogx((cm6{ii,jj,r}.objArray(:,end)), 'r--');    
%             semilogx(rad2deg(cm7{ii,jj,r}.objArray(:,end)), 'k-');
%             semilogx(rad2deg(cm8{ii,jj,r}.objArray(:,end)), 'k--');
            set(gca,'xscale','log');
            grid on;

            title(sprintf('%s / %s', object_names{ii}, graph_names{jj}));
            xlabel('iterations');
            ylabel('global objective');
            legend(method_names, 'Location', 'SouthWest');
            
            saveas(gcf, sprintf('fig/eps/%s_Objective_%s_%s_%02d.eps', ...
                dataset, object_names{idx}, graph_names{jj}, r), 'epsc');
            saveas(gcf, sprintf('fig/png/%s_Objective_%s_%s_%02d.png', ...
                dataset, object_names{idx}, graph_names{jj}, r), 'png');
        end
    end
    if ~vis, close gcf; end
end

% objective (relative)
for r = 1
    figure;
    for jj = range_Network
        for idx = 1:length(range_obj)
            ii=range_obj(idx);
            %subplot(1,length(range_obj),(jj-1)*length(range_obj)+idx);
            clf;
            hold on;
            loglog( ...
                abs(abs(cm1{ii,jj,r}.objArray(2:end,end) - cm1{ii,jj,r}.objArray(1:end-1,end))./cm1{ii,jj,r}.objArray(1:end-1,end)), ...
                'g');
            loglog( ...
                abs(abs(cm2{ii,jj,r}.objArray(2:end,end) - cm2{ii,jj,r}.objArray(1:end-1,end))./cm2{ii,jj,r}.objArray(1:end-1,end)), ...
                'b');
            loglog( ...
                abs(abs(cm3{ii,jj,r}.objArray(2:end,end) - cm3{ii,jj,r}.objArray(1:end-1,end))./cm3{ii,jj,r}.objArray(1:end-1,end)), ...
                'r');
            loglog( ...
                abs(abs(cm4{ii,jj,r}.objArray(2:end,end) - cm4{ii,jj,r}.objArray(1:end-1,end))./cm4{ii,jj,r}.objArray(1:end-1,end)), ...
                'c');
            loglog( ...
                abs(abs(cm5{ii,jj,r}.objArray(2:end,end) - cm5{ii,jj,r}.objArray(1:end-1,end))./cm5{ii,jj,r}.objArray(1:end-1,end)), ...
                'b--');
            loglog( ...
                abs(abs(cm6{ii,jj,r}.objArray(2:end,end) - cm6{ii,jj,r}.objArray(1:end-1,end))./cm6{ii,jj,r}.objArray(1:end-1,end)), ...
                'r--');
%             semilogx(rad2deg(cm7{ii,jj,r}.objArray(:,end)), 'k-');
%             semilogx(rad2deg(cm8{ii,jj,r}.objArray(:,end)), 'k--');
            set(gca,'xscale','log');
            set(gca,'yscale','log');
            grid on;

            title(sprintf('%s / %s', object_names{ii}, graph_names{jj}));
            xlabel('iterations');
            ylabel('relative error of objective');
            legend(method_names, 'Location', 'SouthWest');
            grid on;
            
            saveas(gcf, sprintf('fig/eps/%s_Objective_%s_%s_%02d.eps', ...
                dataset, object_names{idx}, graph_names{jj}, r), 'epsc');
            saveas(gcf, sprintf('fig/png/%s_Objective_%s_%s_%02d.png', ...
                dataset, object_names{idx}, graph_names{jj}, r), 'png');
        end
    end
    if ~vis, close gcf; end
end

% primal norm
for r = 1
    figure;
    for jj = range_Network
        for idx = 1:length(range_obj)
            ii=range_obj(idx);
            %subplot(1,length(range_obj),(jj-1)*length(range_obj)+idx);
            clf;
            hold on;
            semilogx((cm1{ii,jj,r}.rtArray(1:cm1{ii,jj,r}.eITER))./(D*N+D+1), 'g');
            semilogx((cm2{ii,jj,r}.prinormArray(1:cm2{ii,jj,r}.eITER))./(D*N+D+1), 'b');
            semilogx((cm3{ii,jj,r}.prinormArray(1:cm3{ii,jj,r}.eITER))./(D*N+D+1), 'r');    
            semilogx((cm4{ii,jj,r}.prinormArray(1:cm4{ii,jj,r}.eITER))./(D*N+D+1), 'c');
            semilogx((cm5{ii,jj,r}.prinormArray(1:cm5{ii,jj,r}.eITER))./(D*N+D+1), 'b--');
            semilogx((cm6{ii,jj,r}.prinormArray(1:cm6{ii,jj,r}.eITER))./(D*N+D+1), 'r--');     
%             semilogx(rad2deg(cm7{ii,jj,r}.objArray(:,end)), 'k-');
%             semilogx(rad2deg(cm8{ii,jj,r}.objArray(:,end)), 'k--');
            set(gca,'xscale','log');
            set(gca,'yscale','log');
            grid on;

            title(sprintf('%s / %s', object_names{ii}, graph_names{jj}));
            xlabel('iterations');
            ylabel('primal norm');
            legend(method_names, 'Location', 'SouthWest');
            
            saveas(gcf, sprintf('fig/eps/%s_Objective_%s_%s_%02d.eps', ...
                dataset, object_names{idx}, graph_names{jj}, r), 'epsc');
            saveas(gcf, sprintf('fig/png/%s_Objective_%s_%s_%02d.png', ...
                dataset, object_names{idx}, graph_names{jj}, r), 'png');
        end
    end
    if ~vis, close gcf; end
end

% dual norm
for r = 1
    figure;
    for jj = range_Network
        for idx = 1:length(range_obj)
            ii=range_obj(idx);
            %subplot(1,length(range_obj),(jj-1)*length(range_obj)+idx);
            clf;
            hold on;
            semilogx((cm1{ii,jj,r}.stArray(1:cm1{ii,jj,r}.eITER))./(D*N+D+1), 'g');
            semilogx((cm2{ii,jj,r}.dualnormArray(1:cm2{ii,jj,r}.eITER))./(D*N+D+1), 'b');
            semilogx((cm3{ii,jj,r}.dualnormArray(1:cm3{ii,jj,r}.eITER))./(D*N+D+1), 'r');    
            semilogx((cm4{ii,jj,r}.dualnormArray(1:cm4{ii,jj,r}.eITER))./(D*N+D+1), 'c');
            semilogx((cm5{ii,jj,r}.dualnormArray(1:cm5{ii,jj,r}.eITER))./(D*N+D+1), 'b--');
            semilogx((cm6{ii,jj,r}.dualnormArray(1:cm6{ii,jj,r}.eITER))./(D*N+D+1), 'r--');    
%             semilogx(rad2deg(cm7{ii,jj,r}.objArray(:,end)), 'k-');
%             semilogx(rad2deg(cm8{ii,jj,r}.objArray(:,end)), 'k--');
            set(gca,'xscale','log');
            set(gca,'yscale','log');
            grid on;

            title(sprintf('%s / %s', object_names{ii}, graph_names{jj}));
            xlabel('iterations');
            ylabel('dual norm');
            legend(method_names, 'Location', 'SouthWest');
            
            saveas(gcf, sprintf('fig/eps/%s_Objective_%s_%s_%02d.eps', ...
                dataset, object_names{idx}, graph_names{jj}, r), 'epsc');
            saveas(gcf, sprintf('fig/png/%s_Objective_%s_%s_%02d.png', ...
                dataset, object_names{idx}, graph_names{jj}, r), 'png');
        end
    end
    if ~vis, close gcf; end
end

% primal tol
for r = 1
    figure;
    for jj = range_Network
        for idx = 1:length(range_obj)
            ii=range_obj(idx);
            %subplot(1,length(range_obj),(jj-1)*length(range_obj)+idx);
            clf;
            hold on;
            semilogx((cm1{ii,jj,r}.rtArray(1:cm1{ii,jj,r}.eITER))./(D*N+D+1), 'g');
            semilogx((cm2{ii,jj,r}.pritolArray(1:cm2{ii,jj,r}.eITER))./(D*N+D+1), 'b');
            semilogx((cm3{ii,jj,r}.pritolArray(1:cm3{ii,jj,r}.eITER))./(D*N+D+1), 'r');    
            semilogx((cm4{ii,jj,r}.pritolArray(1:cm4{ii,jj,r}.eITER))./(D*N+D+1), 'c');
            semilogx((cm5{ii,jj,r}.pritolArray(1:cm5{ii,jj,r}.eITER))./(D*N+D+1), 'b--');
            semilogx((cm6{ii,jj,r}.pritolArray(1:cm6{ii,jj,r}.eITER))./(D*N+D+1), 'r--');     
%             semilogx(rad2deg(cm7{ii,jj,r}.objArray(:,end)), 'k-');
%             semilogx(rad2deg(cm8{ii,jj,r}.objArray(:,end)), 'k--');
            set(gca,'xscale','log');
            set(gca,'yscale','log');
            grid on;

            title(sprintf('%s / %s', object_names{ii}, graph_names{jj}));
            xlabel('iterations');
            ylabel('primal tolerance');
            legend(method_names, 'Location', 'SouthWest');
            
            saveas(gcf, sprintf('fig/eps/%s_Objective_%s_%s_%02d.eps', ...
                dataset, object_names{idx}, graph_names{jj}, r), 'epsc');
            saveas(gcf, sprintf('fig/png/%s_Objective_%s_%s_%02d.png', ...
                dataset, object_names{idx}, graph_names{jj}, r), 'png');
        end
    end
    if ~vis, close gcf; end
end

% dual tol
for r = 1
    figure;
    for jj = range_Network
        for idx = 1:length(range_obj)
            ii=range_obj(idx);
            %subplot(1,length(range_obj),(jj-1)*length(range_obj)+idx);
            clf;
            hold on;
            semilogx((cm1{ii,jj,r}.stArray(1:cm1{ii,jj,r}.eITER))./(D*N+D+1), 'g');
            semilogx((cm2{ii,jj,r}.dualtolArray(1:cm2{ii,jj,r}.eITER))./(D*N+D+1), 'b');
            semilogx((cm3{ii,jj,r}.dualtolArray(1:cm3{ii,jj,r}.eITER))./(D*N+D+1), 'r');    
            semilogx((cm4{ii,jj,r}.dualtolArray(1:cm4{ii,jj,r}.eITER))./(D*N+D+1), 'c');
            semilogx((cm5{ii,jj,r}.dualtolArray(1:cm5{ii,jj,r}.eITER))./(D*N+D+1), 'b--');
            semilogx((cm6{ii,jj,r}.dualtolArray(1:cm6{ii,jj,r}.eITER))./(D*N+D+1), 'r--');    
%             semilogx(rad2deg(cm7{ii,jj,r}.objArray(:,end)), 'k-');
%             semilogx(rad2deg(cm8{ii,jj,r}.objArray(:,end)), 'k--');
            set(gca,'xscale','log');
            set(gca,'yscale','log');
            grid on;

            title(sprintf('%s / %s', object_names{ii}, graph_names{jj}));
            xlabel('iterations');
            ylabel('dual tolerance');
            legend(method_names, 'Location', 'SouthWest');
            
            saveas(gcf, sprintf('fig/eps/%s_Objective_%s_%s_%02d.eps', ...
                dataset, object_names{idx}, graph_names{jj}, r), 'epsc');
            saveas(gcf, sprintf('fig/png/%s_Objective_%s_%s_%02d.png', ...
                dataset, object_names{idx}, graph_names{jj}, r), 'png');
        end
    end
    if ~vis, close gcf; end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute mean and variance
ssa_all = zeros(length(range_obj),6);
for r = 1
    for jj = 4
        for idx = 1 : length(range_obj)
            ii = range_obj(idx);
            ssa_all(idx,1) = cm1{ii,jj,r}.ssaArray(cm1{ii,jj,r}.eITER);
            ssa_all(idx,2) = cm2{ii,jj,r}.ssaArray(cm2{ii,jj,r}.eITER);
            ssa_all(idx,3) = cm3{ii,jj,r}.ssaArray(cm3{ii,jj,r}.eITER);
            ssa_all(idx,4) = cm4{ii,jj,r}.ssaArray(cm4{ii,jj,r}.eITER);
            ssa_all(idx,5) = cm5{ii,jj,r}.ssaArray(cm5{ii,jj,r}.eITER);
            ssa_all(idx,6) = cm6{ii,jj,r}.ssaArray(cm6{ii,jj,r}.eITER);
        end
    end
end
ssa_all_sorted = sort(ssa_all,1,'ascend'); 
boxplot(rad2deg(ssa_all_sorted(1:(135*0.75),:))); 
ylim([0,100]);
set(gca, 'XTickLabel', {'ADMM', 'ADMM-AP', 'ADMM-NAP', 'ADMM-VP', 'ADMM-VP + AP', 'ADMM-VP + NAP'});
set(gca,'XTickLabelRotation',45);
ylabel('Maximum subspace angle (degrees)');

% compute mean and variance
iter_all = zeros(length(range_obj),6);
for r = 1
    for jj = 4
        for idx = 1 : length(range_obj)
            ii = range_obj(idx);
            iter_all(idx,1) = cm1{ii,jj,r}.eITER;
            iter_all(idx,2) = cm2{ii,jj,r}.eITER;
            iter_all(idx,3) = cm3{ii,jj,r}.eITER;
            iter_all(idx,4) = cm4{ii,jj,r}.eITER;
            iter_all(idx,5) = cm5{ii,jj,r}.eITER;
            iter_all(idx,6) = cm6{ii,jj,r}.eITER;
        end
    end
end
iter_all_sorted = sort(iter_all,1,'ascend'); 
tmp = rad2deg(ssa_all) <= 15;
tmp2 = [mean(iter_all(tmp(:,1),1)) ...
    mean(iter_all(tmp(:,2),2)) ...
    mean(iter_all(tmp(:,3),3)) ...
    mean(iter_all(tmp(:,4),4)) ...
    mean(iter_all(tmp(:,5),5)) ...
    mean(iter_all(tmp(:,6),6))];
bar(tmp2); 
set(gca, 'XTickLabel', {'ADMM', 'ADMM-AP', 'ADMM-NAP', 'ADMM-VP', 'ADMM-VP + AP', 'ADMM-VP + NAP'});
set(gca,'XTickLabelRotation',45);
ylabel('Number of iterations');
