%--------------------------------------------------------------------------
% SfM Hopkins155 data demonstration
%
% Implemented/Modified
%  by     Sejong Yoon (sjyoon@cs.rutgers.edu)
%  on     2015.04.03 (last modified on 2015/04/03)
%--------------------------------------------------------------------------

% Clear data
clear; close all;

% Choose random seed: optional setting to reproduce numbers. You should be
% able to obtain consistent result without this although numbers may be
% different from those reported in the paper.
s = RandStream('mt19937ar','Seed',1);
RandStream.setGlobalStream(s);
reset(s,0);

% Create result output directory if needed
if ~exist('results/sfm_hopkins', 'dir');
    mkdir('results/sfm_hopkins');
end

%% Options
objs = 135; % total number of objects
repeats = 5; % repeat this number of individual trials

M = 3; % latent space dimension

THRESHc = 1e-3; % convergence precision
THRESHd = 1e-3; % convergence precision

ETA = 10; % learning parameter
NV = 5; % number of nodes

% graph topology (1 = compete, 2 = ring)
E = get_adj_graph(NV); E = E{2}; 

% Frequency for intermediate object function value output
objfreq_c = 0;
objfreq_d = 0;

%% For each video, we run all methods and save angles...
SAs = cell(6, 1);
for j = 1 : 6
    SAs{j} = zeros(objs, repeats);
end

for idm = 1 : objs
    tmp = load(sprintf('data/hopkins/%03d.mat', idm));
    
    fprintf('[[ Current Object: %d ]]\n', idm);
    mm_mat = tmp.measurement_matrix;

    % get dimension of measurement matrix
    [D, N] = size(mm_mat);

    % Measurement matrix should be in the specific form
    if mod(D, 2) ~= 0
        error('Measurement matrix should be (2 x #frames)x(# points) form!');
    end
    
    % node assignment of samples(=frames)
    V = get_sample_assign(NV, D, 'IsSfM', true);

    % translate origin of measurement matrix to zero
    centroid = mean(mm_mat, 2);
    mm_trans = mm_mat - repmat(centroid, [1, N]);
    X = mm_trans';

    % We consider centralized SVD SfM result as ground truth reference
    [sfm_M, sfm_S, sfm_b, ~, ~, ~] = sfm_affine(mm_mat);
    GT = sfm_S'; % We consider SVD SfM as ground truth
    
    for idr = 1 : repeats % independent runs (initial value)
        % initialize
        m_init_c = get_init_value(X, M, 'ModelType', 'cppca');
        m_init_d = get_init_value(X, M, 'ModelType', 'sfm_d', ...
            'NumberOfNodes', NV, 'SampleAssignVec', V);  

        % run
%         [cm1, cm2] = expr_run_cppca(X, M, THRESHc, objfreq_c, m_init_c, true);
%         [cm3, cm4] = expr_run_dppca(X, M, V, E, ETA, THRESHd, objfreq_d, m_init_d, true);
%         [cm5, cm6] = expr_run_vbpca(X, M);
        load(sprintf('results/sfm_hopkins/%03d_%03d.mat', idm, idr));
        
        % Centralized settting result
        angle1 = subspace(GT, cm1.W);
%         angle2 = subspace(GT, cm2.mW);

%         angle5 = subspace(GT, cm5.W);
%         angle6 = subspace(GT, cm6.mW);
        
        % Distributed settting result
        angle3 = subspace(GT, cm3.W{1});
        angle4 = subspace(GT, cm4.W{1});
        for idj = 2:NV
            if angle3 < subspace(GT, cm3.W{idj})
                angle3 = subspace(GT, cm3.W{idj});
            end
            if angle4 < subspace(GT, cm4.W{idj})
                angle4 = subspace(GT, cm4.W{idj});
            end
        end
        
        % iterations
%         angle1 = cm1.eITER;
%         angle3 = cm3.eITER;
%         angle4 = cm4.eITER;

        % time
        angle1 = cm1.eTIME;
        angle3 = cm3.eTIME;
        angle4 = cm4.eTIME;        

        angles = [angle1 0 angle3 angle4 0 0];
%         angles = [angle1 angle2 angle3 angle4 angle5 angle6];
        
        % Print Result
        fprintf('*** Result ***\n');
        fprintf(' ### Subspace angle (S/A)                 :1e-123456789    15\n');
        
        fprintf('(1) GT and (2) PPCA                   : %.15f (Init: %02d / %f)\n', ...
            angle1, idr, m_init_c.VAR);
%         fprintf('(1) GT and (2) BPCA                   : %.15f (Init: %02d / %f)\n', ...
%             angle2, idr, m_init_c.VAR);

        fprintf('(1) GT and (2) D-PPCA    (Nodes:%d)   : %.15f (Init: %02d / %f)\n', ...
            NV, angle3, idr, m_init_d.VAR{1});
        fprintf('(1) GT and (2) D-PPCA-ANT(Nodes:%d)   : %.15f (Init: %02d / %f)\n', ...
            NV, angle4, idr, m_init_d.VAR{1});

%         fprintf('(1) GT and (2) PPCA (IR)              : %.15f \n', ...
%             angle5, idr);
%         fprintf('(1) GT and (2) BPCA (IR)              : %.15f \n', ...
%             angle6, idr);
        
        parfor_var_save9( ...
            sprintf('results/sfm_hopkins/%03d_%03d.mat', idm, idr), ...
            cm1, cm2, cm3, cm4, cm5, cm6, m_init_c, m_init_d, angles);
    end % parfor: repeats independent runs (initial value)

    %----------------------------------------------------------
    % Store Result
    for idr = 1 : repeats
        tmp = load(sprintf('results/sfm_hopkins/%03d_%03d.mat', idm, idr));
        for j = 1 : 6
            SAs{j}(idm, idr) = tmp.angles(j);
        end
    end

    %--------------------------------------------------------------
    % Save the model for future usage
    ms_path = sprintf('results/sfm_hopkins/%03d_svd.mat', idm);
    save(ms_path, 'sfm_M', 'sfm_S', 'sfm_b');
end

clearvars -except SAs;

%% Output overall result
% save('results/sfm_hopkins/all_results.mat','SAs');
% 
% load('results/sfm_hopkins/all_results.mat');

fprintf('\n\n');
method = {'PPCA', 'BPCA', 'D-PPCA', 'D-PPCA-ANT', 'PPCA (IR)', 'BPCA (IR)'};
for j = [1,3,4]%1 : 6
    % mean over independent runs
    angle_mean = mean(mean(SAs{j}));% * 180 / pi;
    angle_std = std(std(SAs{j}, 1, 2), 1, 1);% * 180 / pi;

    fprintf('### S/A b/w GT and %s): %f (degrees) +- %f (degree)\n', ...
        method{j}, angle_mean, angle_std);
end
