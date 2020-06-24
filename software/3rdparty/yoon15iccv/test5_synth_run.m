%--------------------------------------------------------------------------
% Synthetic data demonstration
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

% Path
addpath(genpath('PCAMV'));

% create data and result directory if needed
if ~exist('data', 'dir'), mkdir('data'); end
if ~exist('data/synth', 'dir'), mkdir('data/synth'); end
if ~exist('results', 'dir'), mkdir('data'); end
if ~exist('results/synth', 'dir'), mkdir('results/synth'); end

%% ------------------------------------------------------------------------
% Experiment setting

% V: Node assignment of samples
Varr = {1, 2, 5, 8, 10};

% ETA: Learning rate
ETAarr = {8, 10, 12, 16};

% Missing rate and repetition
RATEarr = {0.01, 0.05, 0.1, 0.2, 0.3};
REPEATS = 10;

% Number of nodes and ETA value used for missing value experiments; these
% are indices to Varr and ETAarr. e.g., MISS_V = 3 means V = 5.
MISS_V = 3;
MISS_ETA = 2;
MISS_E = 2;

% Miscellaneous settings
THRESHc = 1e-2;     % threshold for centralized methods (PPCA, BPCA)
THRESHd = 1e-2;     % threshold for distributed methods (D-PPCA, D-BPCA)
objfreq_c = 10;     % show objective per this iterations (PPCA, BPCA)
objfreq_d = 10;     % show objective per this iterations (D-PPCA, D-BPCA)

% Save options for record
save('data/synth/options.mat', ...
    'Varr', 'ETAarr', 'RATEarr', 'REPEATS', 'THRESHc', 'THRESHd', ...
    'objfreq_c', 'objfreq_d', 'MISS_V', 'MISS_ETA', 'MISS_E');

%% ------------------------------------------------------------------------
% Generate data

% 20 dimensional, 500 samples searching, 5 dimensional subspace
VarW = rand(1);
VarMU = rand(1);
VarX = rand(1);

mW = 10;
mMU = 5;
mX = 0;

N = 500;
D = 20;
M = 5;

% Z  (dim=M)   comes from N(0,I) 
Z  = mvnrnd(zeros(N,M), eye(M)); % each row is a sample
% W  (dim=D,M) comes from N(mW,vW)
W = orth(mvnrnd(ones(D,M)*mW, eye(M)*VarW)) * diag(M:-1:1);
% MU (dim=D)   comes from N(mMU,vMU)
MU = repmat(mvnrnd(ones(1,D)*mMU, eye(D)*VarMU), [N,1]);
% E  (dim=D)   comes from N(mX,vX)
%E = mvnrnd(ones(N,D)*mX, eye(D)*VarX);
E = ones(N,D)*mX + VarX*randn(N, D);
% model
X = W * Z' + MU' + E';

% save generated data for record
save('data/synth/data.mat', ...
    'VarX', 'VarMU', 'VarX', 'mW', 'mMU', 'mX', ...
    'N', 'D', 'M', 'Z', 'W', 'MU', 'E', 'X');

% Missing value data
parfor r = 1 : length(RATEarr)
    MissRate = RATEarr{r};
    
    for RunIdx = 1 : REPEATS
        % generate MAR corrupted data
        Xmar = X;
        seq = randperm(D * N);
        seq = seq(1:floor(D * N * MissRate / 100));
        Xmar(seq) = NaN;
        
        % generate MNAR corrupted data
        Xmnar = X;
        seq = sfm_gen_missing_mnar(X, MissRate);
        Xmnar(seq == 0) = NaN;
        
        parfor_var_save3(sprintf('data/synth/data_R%d_i%02d.mat', r, RunIdx), ...
            Xmar, Xmnar, MissRate);
    end
end

%% ------------------------------------------------------------------------
% Initializations

% cent: initializes with a random projection, zero mean and a small variance
m_init_c = get_init_value(X, M, 'ModelType', 'cppca');

% dist: initializes with local PPCA (SVD)
m_init_d = cell(length(Varr), 1);
for idk = 1 : length(Varr)
    % Node assignment to each sample
    NV = Varr{idk};
    V = get_sample_assign(NV, N);

    m_init_d{idk} = get_init_value(X, M, 'ModelType', 'dppca', ...
        'NumberOfNodes', NV, 'SampleAssignVec', V);    
    
    fprintf('Generated initialization for NV = %d case.\n', NV);
end

% dist: initializes with local PPCA with missing values (EM)
NV = Varr{MISS_V}; 
V = get_sample_assign(NV, N);

parfor k = 1 : length(RATEarr)
    for l = 1 : REPEATS
        tmp = parfor_var_load(sprintf('data/synth/data_R%d_i%02d.mat', k, l));

        m_init_d_mar = get_init_value(tmp.Xmar, M, 'ModelType', 'dppca', ...
            'NumberOfNodes', NV, 'SampleAssignVec', V);

        m_init_d_mnar = get_init_value(tmp.Xmnar, M, 'ModelType', 'dppca', ...
            'NumberOfNodes', NV, 'SampleAssignVec', V);
        
        fprintf('Generated initialization for noise %f (%d/%d).\n', ...
            RATEarr{k}, l, REPEATS);

        parfor_var_save2(sprintf('data/synth/init_R%d_i%02d.mat', k, l), ...
            m_init_d_mar, m_init_d_mnar);
    end
end

% save initializations for record
save('data/synth/init.mat', 'm_init_c', 'm_init_d');

%% ------------------------------------------------------------------------
% Experiments 1

% Various number of nodes, network topology and ETA
% NOTE: We used parallel computation toolbox of MATLAB here but one can
%  simply change [parfor] with [for] to run the code without the toolbox.
parfor idk = 1:length(Varr)
    % Node assignment to each sample
    NV = Varr{idk};
    V = get_sample_assign(NV, N);
    
    % Adjacency graph topology (1:complete, 2:ring, 3:star, 4:chain)
    GRAPHS = get_adj_graph(NV);
    
    for idx = 1:length(GRAPHS)
        % choose one graph type
        E = GRAPHS{idx};

        for idy = 1:length(ETAarr)
            % choose ETA
            ETA = ETAarr{idy};

            [cm1, cm2] = expr_run_cppca(X, M, THRESHc, objfreq_c, m_init_c, false);
            [cm3, cm4] = expr_run_dppca(X, M, V, E, ETA, THRESHd, objfreq_d, m_init_d{idk}, false);
            [cm5, cm6] = expr_run_vbpca(X, M);
            
            savepath = sprintf('results/synth/N%02d_G%03d_E%02d.mat', NV, idx, ETA);
            parfor_var_save6(savepath, cm1, cm2, cm3, cm4, cm5, cm6);
        end
    end
end

%% ------------------------------------------------------------------------
% Prepare missing value experiments
NV = Varr{MISS_V}; 
V = get_sample_assign(NV, N);
ETA = ETAarr{MISS_ETA};
GRAPHS = get_adj_graph(NV); 
E = GRAPHS{MISS_E};

%% ------------------------------------------------------------------------
% Experiments 2 (MAR)
parfor k = 1 : length(RATEarr)
    for l = 1 : REPEATS
        tmp1 = parfor_var_load(sprintf('data/synth/data_R%d_i%02d.mat', k, l));
        tmp2 = parfor_var_load(sprintf('data/synth/init_R%d_i%02d.mat', k, l));

        [cm1, cm2] = expr_run_cppca(tmp1.Xmar, M, THRESHc, objfreq_c, m_init_c, false);
        [cm3, cm4] = expr_run_dppca(tmp1.Xmar, M, V, E, ETA, THRESHd, objfreq_d, tmp2.m_init_d_mar, false);
        [cm5, cm6] = expr_run_vbpca(tmp1.Xmar, M);

        savepath = sprintf('results/synth/R%d_i%02d_mar.mat', k, l);
        parfor_var_save6(savepath, cm1, cm2, cm3, cm4, cm5, cm6);
    end
end

%% ------------------------------------------------------------------------
% Experiments 3 (MNAR)
parfor k = 1 : length(RATEarr)
    for l = 1 : REPEATS
        tmp1 = parfor_var_load(sprintf('data/synth/data_R%d_i%02d.mat', k, l));
        tmp2 = parfor_var_load(sprintf('data/synth/init_R%d_i%02d.mat', k, l));

        [cm1, cm2] = expr_run_cppca(tmp1.Xmnar, M, THRESHc, objfreq_c, m_init_c, false);
        [cm3, cm4] = expr_run_dppca(tmp1.Xmnar, M, V, E, ETA, THRESHd, objfreq_d, tmp2.m_init_d_mnar, false);
        [cm5, cm6] = expr_run_vbpca(tmp1.Xmar, M);

        savepath = sprintf('results/synth/R%d_i%02d_mnar.mat', k, l);
        parfor_var_save6(savepath, cm1, cm2, cm3, cm4, cm5, cm6);
    end
end

%% ------------------------------------------------------------------------
% Plot result
test5_synth_show;
