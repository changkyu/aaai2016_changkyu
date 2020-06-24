%% ------------------------------------------------------------------------
% Synthetic data demonstration (w/ or w/o missing values)

%% ------------------------------------------------------------------------
% Generate data

clear; close all;

% Choose random seed: optional setting to reproduce numbers.
s = RandStream('mt19937ar','Seed',1);
RandStream.setGlobalStream(s);
reset(s,0);

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

%% ------------------------------------------------------------------------
% Experiment setting

% ETA: Learning rate
ETA = 10;

% Node assignment to each sample
NV = 2;
V = get_sample_assign(NV, N);

% E: Adjacency graph topology (1:complete, 2:ring, 3:star, 4:chain)
GRAPHS = get_adj_graph(NV);
E = GRAPHS{1};

% Miscellaneous settings
THRESHc = 1e-3;     % threshold for centralized methods (PPCA, VBPCA)
THRESHd = 1e-3;     % threshold for distributed methods (D-PPCA, D-BPCA)
objfreq_c = 1;      % show objective per this iterations (PPCA, VBPCA)
objfreq_d = 1;      % show objective per this iterations (D-PPCA, D-BPCA)

% Missing rate of data
MissRate = 0;       % missing rate in percentage
if MissRate > 0
    seq = randperm(D * N);
    seq = seq(1:floor(D * N * MissRate / 100));
    X(seq) = NaN;
end

%% ------------------------------------------------------------------------
% Run an experiment

% initializes with a random projection, zero mean and a small variance
m_init_c = get_init_value(X, M, 'ModelType', 'cppca');

[cm1, cm2] = expr_run_cppca(X, M, THRESHc, objfreq_c, m_init_c, false);

% initializes with local PPCA (either SVD or EM)
NV = max(V); % number of nodes
m_init_d = get_init_value(X, M, 'ModelType', 'dppca', ...
    'NumberOfNodes', NV, 'SampleAssignVec', V);

[cm3, cm4] = expr_run_dppca(X, M, V, E, ETA, THRESHd, objfreq_d, m_init_d, false);

%% ------------------------------------------------------------------------
% Report
% res = expr_calc_comparison(cm1, cm2, cm3, cm4, ...
%     X, W, Z, MU, mX, mW, mMU, VarX, VarW, VarMU);

