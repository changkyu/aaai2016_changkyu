%--------------------------------------------------------------------------
% Synthetic data demonstration
%
% Implemented/Modified
%  by     Sejong Yoon (sjyoon@cs.rutgers.edu)
%  on     2012.02.05 (last modified on 2015/02/08)
%--------------------------------------------------------------------------

% Clear data
clear;

% Frequency for intermediate object function value output
objfreq_c = 100;
objfreq_d1 = 100;
objfreq_dj = 100;

% Choose random seed: optional setting to reproduce numbers. You should be
% able to obtain consistent result without this although numbers may be
% different from those reported in the paper. In some cases, the D-PPCA
% objective function value may fluctuate near the stationary point and it
% may take a while to stabilize. However, they will converge in the end.
s = RandStream('mt19937ar','Seed',1);
RandStream.setGlobalStream(s);
reset(s,0);

% 50 dimensional, 100 samples searching, 5 dimensional subspace
N = 100; D = 50; M = 5; VAR = 1;
% Z (dim=M) comes from N(0,I)
Z = randn(s,N,M);
% W (dim=M) comes from N(0,I)
% NOTE: W can be an arbitrary matrix, i.e. W = rand(s,D,M);
W = randn(s,D,M);
% E (dim=D) comes from N(0,VAR*I)
E = (sqrt(VAR) .* randn(s,N,D))';
% Our PPCA model is X = W * Z + E
measurement_matrix = W * Z' + E;
datastr = '';

% get dimension of measurement matrix
[D, N] = size(measurement_matrix);

mm_trans = measurement_matrix;

% create result directory if needed
if ~exist('result/synth_dppca', 'dir')
    mkdir('result/synth_dppca');
end

%% ------------------------------------------------------------------------
disp('*** Centralized Setting ***');

M = 5; % latent space dimension
THRESHc = 1e-5; % convergence precision (default: 1e-5)

disp('* PPCA (SVD) *');
cm1 = cppca(mm_trans, M);

% initializes with a random projection, zero mean and a small variance
m_init1 = get_init_value(mm_trans, M, ...
    'ModelType', 'cppca', 'VarianceFactor', 1);

disp('* PPCA (EM) *');
cm2 = cppca_em(mm_trans, M, ...
    'Threshold', THRESHc, 'InitModel', m_init1, 'ShowObjPer', objfreq_c);

%% ------------------------------------------------------------------------
% Experiments
disp('*** Distributed Setting ***');

THRESHd = 1e-5; % convergence precision (default: 1e-5)

% V: Node assignment of samples
Varr = {1, 2, 5, 8, 10};

% ETA: Learning rate
ETAarr = {8, 10, 12, 16};

% Various number of nodes, network topology and ETA
% NOTE: We used parallel computation toolbox of MATLAB here but one can
%  simply change [parfor] with [for] to run the code without the toolbox.
parfor idk = 1:length(Varr)
    % Node assignment to each sample
    NV = Varr{idk};
    Vp = get_sample_assign(NV, N);
    Vc = get_sample_assign(NV, D);

    m_init3 = get_init_value(mm_trans, M, 'ModelType', 'dppca_dup', ...
        'NumberOfNodes', NV, 'VarianceFactor', 10, 'Init', m_init1);
    
    % E: Network topology
    Earr = get_adj_graph(NV);

    for idx = 1:length(Earr)
        E = Earr{idx};

        for idy = 1:length(ETAarr)
            ETA = ETAarr{idy};

            cm3 = dppca(mm_trans,  M, Vp, E, 'InitModel', m_init3, ...
               'Eta', ETA, 'Threshold', THRESHd, 'ShowObjPer', objfreq_dj);
            
            savepath = sprintf('result/synth_dppca/dppca_N%02d_G%03d_E%02d.mat', NV, idx, ETA);
            save_var_parfor1(savepath, cm3);
        end
    end
end

%% ------------------------------------------------------------------------
% Plot result
show_result_1_synth_dppca;
