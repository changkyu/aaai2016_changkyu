%--------------------------------------------------------------------------
% Synthetic data demonstration (Missing Values)
%
% Implemented/Modified
%  by     Sejong Yoon (sjyoon@cs.rutgers.edu)
%  on     2012.06.01 (last modified on 2015/02/08)
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
N = 100; D = 50; M = 5; VAR = 0.2;
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
if ~exist('result/synth_dppca_m', 'dir')
    mkdir('result/synth_dppca_m');
end

% Set thresholds
THRESHc = 1e-5; % convergence precision (default: 1e-5)
THRESHd = 1e-5; % convergence precision (default: 1e-5)

%% MAR --------------------------------------------------------------------
disp('*** Different missing value rate (MAR) ***');

RATEArr = {0.01, 0.05, 0.1, 0.2, 0.3};

parfor idk = 1:length(RATEArr)
    RATE = RATEArr{idk};

    % We simply remove values randomly (1 = missing)
    MissIDX = zeros(D, N);
    MissIDX(randi(D*N, 1, ceil(D*N*RATE))) = 1;

    savepath = sprintf('result/synth_dppca_m/data_R%d.mat', idk);
    save_var_parfor2(savepath, mm_trans, MissIDX);

    for idr = 1:5%1:20 % 20 is too many
        % Run centralized
        m_init1 = get_init_value(mm_trans, M, ...
            'ModelType', 'cppca', 'VarianceFactor', 1);

        fprintf('* PPCA (EM) / %f *\n', RATE);
        cm1 = cppca_em(mm_trans, M, 'Threshold', THRESHc, ...
            'InitModel', m_init1, 'ShowObjPer', objfreq_c, ...
            'MissIDX', MissIDX);

        savepath = sprintf('result/synth_dppca_m/cppca_R%d_i%02d.mat', idk, idr);
        save_var_parfor2(savepath, cm1, MissIDX);

        % Run distributed
        ETA = 10;
        NV = 5;
        Earr = get_adj_graph(NV);
        E = Earr{2};
        Vp = get_sample_assign(NV, N);

        % To remove initialization effect, we duplicate centralized
        % initialization parameters; Note that our interest here is dealing
        % with missing values, not initialization. Moreover, it is
        % possible for all nodes to have the same initialization values.
        m_init2 = get_init_value(mm_trans, M, ...
            'ModelType', 'dppca_dup', 'NumberOfNodes', NV, ...
            'VarianceFactor', 10, 'Init', m_init1);

        cm2 = dppca(mm_trans, M, Vp, E, 'InitModel', m_init2, ...
               'Eta', ETA, 'Threshold', THRESHd, 'ShowObjPer', objfreq_dj, ...
               'MissIDX', MissIDX);

        savepath = sprintf('result/synth_dppca_m/dppca_R%d_i%02d.mat', idk, idr);
        save_var_parfor1(savepath, cm2);
    end
end

%% MNAR -------------------------------------------------------------------
disp('*** Different missing value rate (MNAR) ***');

RATEArr = {0.01, 0.05, 0.1, 0.2, 0.3};

parfor idk = 1:length(RATEArr)
    RATE = RATEArr{idk};

    % We get band matrix: correlated missing values at non-diagonal
    MissIDX = gen_sfm_missing_mnar(mm_trans, RATE);

    savepath = sprintf('result/synth_dppca_m/data_NR%d.mat', idk);
    save_var_parfor2(savepath, mm_trans, MissIDX);

    for idr = 1:5%1:20 % 20 is too many
        % Run centralized
        m_init1 = get_init_value(mm_trans, M, ...
            'ModelType', 'cppca', 'VarianceFactor', 1);

        fprintf('* PPCA (EM) / %f *\n', RATE);
        cm1 = cppca_em(mm_trans, M, 'Threshold', THRESHc, ...
            'InitModel', m_init1, 'ShowObjPer', objfreq_c, ...
            'MissIDX', MissIDX);

        savepath = sprintf('result/synth_dppca_m/cppca_NR%d_i%02d.mat', idk, idr);
        save_var_parfor2(savepath, cm1, MissIDX);

        % Run distributed
        ETA = 10;
        NV = 5;
        Earr = get_adj_graph(NV);
        E = Earr{2};
        Vp = get_sample_assign(NV, N);

        % To remove initialization effect, we duplicate centralized
        % initialization parameters; Note that our interest here is dealing
        % with missing values, not initialization. Moreover, it is
        % possible for all nodes to have the same initialization values.
        m_init2 = get_init_value(mm_trans, M, ...
            'ModelType', 'dppca_dup', 'NumberOfNodes', NV, ...
            'VarianceFactor', 10, 'Init', m_init1);
 
        cm2 = dppca(mm_trans, M, Vp, E, 'InitModel', m_init2, ...
               'Eta', ETA, 'Threshold', THRESHd, 'ShowObjPer', objfreq_dj, ...
               'MissIDX', MissIDX);

        savepath = sprintf('result/synth_dppca_m/dppca_NR%d_i%02d.mat', idk, idr);
        save_var_parfor1(savepath, cm2);
    end
end

%% ------------------------------------------------------------------------
% Plot result
show_result_2_synth_dppca_m
