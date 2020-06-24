%--------------------------------------------------------------------------
% SfM Synthetic data demonstration (MAR and MNAR)
%
% Implemented/Modified
%  by     Sejong Yoon (sjyoon@cs.rutgers.edu)
%  on     2012.02.05 (last modified on 2015/04/02)
%--------------------------------------------------------------------------

% Clear data
clear; close all;

% Choose random seed: optional setting to reproduce numbers. You should be
% able to obtain consistent result without this although numbers may be
% different from those reported in the paper.
s = RandStream('mt19937ar','Seed',1);
RandStream.setGlobalStream(s);
reset(s,0);

model = 'cube';

missings = 10; % Number of missing values samples (MAR)
repeats = 20; % Number of independent runs of different initial values
miss_rate = 0.2; % Fixed missing rate

% Results are stored here with each row with different missing value set
SAs = cell(6, 1);
for j = 1 : 6
    SAs{j} = zeros(missings+1, repeats);
end

%% ------------------------------------------------------------------------
% Data

% Our data is fixed this time. We have to remove the first point from the
% cube data because it is occluded over all frames, thus never observed.
% This will not harm the validity of our result since this will only
% increase difficulties for MNAR settings. We use all points for MAR since
% all points have same probability of missing.
load(sprintf('data/%s/data_missing.mat', model));

% get dimension of measurement matrix
[D, N] = size(measurement_matrix);

% Measurement matrix should be in the specific form
if mod(D, 2) ~= 0
    error('Measurement matrix should be (2 x #frames)x(# points) form!');
end

% Prepare missing value matrices
%   First MissIDX = MNAR; 
%   Other missing value sets (indices) = MAR;
MissIDXs = cell(missings+1,1);
MissIDXs{1} = MissIDX; % this is MNAR 
for idx = 2:(missings+1)
    MissIDXs{idx} = sfm_gen_missing_mar(measurement_matrix, miss_rate);
end

% Create result output directory if needed
if ~exist(sprintf('results/sfm_%s_m', model), 'dir');
    mkdir(sprintf('results/sfm_%s_m', model));
end

% For missing MAR cases plus one MNAR case...
for idm = 1 : length(MissIDXs)
    % Pick current missing value set (now 1 = available, 0 = missing)
    MissIDX = logical(MissIDXs{idm});
    
    % translate origin of measurement matrix to zero
    centroid = sum(measurement_matrix.*MissIDX, 2) ./ sum(MissIDX, 2);
    if idm == 1
        mm_trans = measurement_matrix - repmat(centroid, [1, N]);
    else
        mm_trans = measurement_matrix - repmat(centroid, [1, N]);
    end
    
    X = mm_trans;
    X(~MissIDX) = NaN;
    X = X';
    
    %% Options

    M = 3; % latent space dimension

    THRESHc = 1e-3; % convergence precision
    THRESHd = 1e-3; % convergence precision
    
    ETA = 10; % learning parameter
    NV = 5; % number of nodes

    % Frequency for intermediate object function value output
    objfreq_c = 0;
    objfreq_d = 0;
    
    %% --------------------------------------------------------------------
    % PPCA, BPCA, D-PPCA, D-BPCA

    % node assignment of samples(=frames)
    V = get_sample_assign(NV, D, 'IsSfM', true);
    % graph topology (1 = compete, 2 = ring)
    E = get_adj_graph(NV); E = E{2}; 
    
    parfor idr = 1:repeats % independent runs (initial value)
        % initialize
        m_init_c = get_init_value(X, M, 'ModelType', 'cppca');
        m_init_d = get_init_value(X, M, 'ModelType', 'sfm_d', ...
            'NumberOfNodes', NV, 'SampleAssignVec', V);  

        % run
        [cm1, cm2] = expr_run_cppca(X, M, THRESHc, objfreq_c, m_init_c, true);
        [cm3, cm4] = expr_run_dppca(X, M, V, E, ETA, THRESHd, objfreq_d, m_init_d, true);
        [cm5, cm6] = expr_run_vbpca(X, M);

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

        angles = [angle1 0 angle3 angle4 0 0];
%         angles = [angle1 angle2 angle3 angle4 angle5 angle6];

        parfor_var_save9( ...
            sprintf('results/sfm_%s_m/%03d_%03d.mat', model, idm, idr), ...
            cm1, cm2, cm3, cm4, cm5, cm6, m_init_c, m_init_d, angles);
    end

    %% Store Result
    for idr = 1:repeats
        tmp = load(sprintf('results/sfm_%s_m/%03d_%03d.mat', model, idm, idr));
        for j = 1 : 6
            SAs{j}(idm, idr) = tmp.angles(j);
        end
    end
end

%% Output overall result
save(sprintf('results/sfm_%s_m/%s_all.mat', model, model), 'SAs');

load(sprintf('results/sfm_%s_m/%s_all.mat', model, model));

fprintf('\n\n');
method = {'PPCA', 'BPCA', 'D-PPCA', 'D-PPCA-ANT', 'PPCA (IR)', 'BPCA (IR)'};
for j = [1,3,4]%1 : 6
    mean_dist = mean(SAs{j}, 2);

    angle_mar = mean(mean_dist(2:end)) * 180 / pi;
    angle_mnar = mean_dist(1) * 180 / pi;

    fprintf('### Result MAR  (s/a b/w g/t and %s): %f (degrees)\n', method{j}, angle_mar);
    fprintf('### Result MNAR (s/a b/w g/t and %s): %f (degrees)\n', method{j}, angle_mnar);
end
