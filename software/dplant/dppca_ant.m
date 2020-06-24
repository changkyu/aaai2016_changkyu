function model = dppca_ant( X, M, V, E, varargin )
% DPPCA_ANT   Distributed Probablistic PCA with Adaptive Penalty
% 
% Description
%  Solve probabilistic PCA problem in a distributed way. The network has 
%  max(V) nodes. We assume the network is connected. This function only 
%  simulates parameter broadcasts. Local computation is done by dppca_local 
%  function. NaN elements in X are considered as missing values.
%
% Input
%  X     : D x N matrix for full data from all nodes (N=N1+N2+...+NJ)
%  M     : Scalar of projection dimension
%  V     : N x 1 vector for each observation's source (node affiliation)
%  E     : J x J adjacency matrix where J = max(V)
%  [Optional Parameters]
%  InitModel  : D-PPCA model to set initial parameter (Def: random)
%  Threshold  : Scalar convergence criterion (objective) (Def: 1e-5)
%  ThresholdA : Scalar convergence criterion (absolute)  (Def: 1e-4)
%  ThresholdR : Scalar convergence criterion (relative)  (Def: 1e-2)
%  ShowObjPer : If > 0, print out objective every specified iteration.
%               If 0, nothing will be printed. (Def: 1)
%  MaxIter    : Maximum iterations (Def: 10000)
%  ZeroMean   : True if we enforce the mean to be zero. (Def: false)
%  Eta        : Scalar of learning rate (Def: 10)
%  VaryETA    : Varying Eta value as in (He, et al. 2000)
%  VaryETA_mu : Parameter mu  for (He, et al. 2000) (Def: 0.1) (0 < mu < 1)
%  VaryETA_tau: Parameter tau for (He, et al. 2000) (Def: 1)   (0 < tau)
%  VaryETA_max: Max iterations    (He, et al. 2000) (Def: 50)
%  W_GT       : Ground truth projection matrix if you want iter-to-iter
%               subspace angle history (Def: [])
%
% Output
%  model = structure(W, MU, VAR, ...);
%  W        : J cells; D x M projection matrices for J nodes
%  MU       : J cells; D x 1 vector sample means for J nodes
%  VAR      : J cells; Scalar estimated variances for J nodes
%  EZ       : J cells; M x N matrix, mean of N latent vectors
%  EZZt     : J cells; M x M x N cube, covariance of N latent vectors
%  eITER    : Iterations took
%  eTIME    : Elapsed time
%  objArray : Objective function value change over iterations
%  LAMBDAi  : J cells; D x M matrix Lagrange multipliers
%  GAMMAi   : J cells; D x 1 vector Lagrange multipliers
%  BETAi    : J cells; Scalar       Lagrange multipliers
%
% Implemented
%  by     Changkyu Song (changkyu.song@rutgers.edu)
%         Sejong Yoon   (sjyoon@cs.rutgers.edu)
%  on     2011.10.07 (last modified on 2015/04/12)
%
% References
%  [1] M.E. Tipping and C.M. Bishop, Probablistic principal component 
%      analysis, J. Royal Statistical Society B 21(3), pp. 611-622, 1999.
%  [2] Probablistic Modeling Toolkit 3, pmtk3.googlecode.com
%  [3] P.A. Forero, A. Cano and G.B. Giannakis, Distributed clustering
%      using wireless sensor networks, IEEE J. Selected Topics in Signal 
%      Processing 5(4), August 2011.
%  [4] S. Yoon and V. Pavlovic, Distributed Probabilistic Learning for 
%      Camera Networks with Missing Data, NIPS 25, 2012.
%  [5] B.S. He and H. Yang and S.L. Wang, Alternating Direction Method with 
%      Self-Adaptive Penalty Parameters for Monotone Variational 
%      Inequalities, J. of Opt. Th. and App. 106(2), pp. 337-356, 2000.
%  [6] T. Goldstein, B. O'Donoghue, S. Setzer and R. Baraniuk, Fast 
%      Alternating Direction Optimization Methods, SIAM J. Imaging Science 
%      7(3) pp. 1588--1623, 2014.

% Check required arguments
assert(nargin >= 4, 'Please specify at least X, M, V and E.');

% D dimensions x N samples
[D, N] = size(X);

% J = number of nodes
J = max(V);

% Check graph is valid
[r,c] = size(E);
assert(r == c, 'Adjacency matrix is not square!');
assert(sum(sum(abs(E' - E))) == 0, 'Graph should be indirectional!');
if J > 1
    assert(r == J, 'Adjacency matrix size does not match number of nodes!');
end

%--------------------------------------------------------------------------
% Parse optional parameters
p = inputParser;
p.StructExpand = false;

W = cell(J,1);
MU = cell(J,1);
VAR = cell(J,1);
for j = 1 : J
    W{j} = orth(randn(D, M));
    MU{j} = zeros(D, 1);
    VAR{j} = 1;
end
defaultMODEL = structure(W, MU, VAR);
defaultTHRESH = 1e-5;
defaultTHRESHA = 1e-4;
defaultTHRESHR = 1e-2;
defaultITER = 1;
defaultMaxIter = 10000;
defaultZeroMean = false;
defaultETA = 10;
defaultVaryETA = false;
defaultVaryETA_mu = 0.1;
defaultVaryETA_tau = 1;
defaultVaryETA_max = 50;
defaultUseResidual = false;
defaultW_GT = [];
defaultFwRes = false;
defaultFwResEta = 0.999;

defaultModelType = 'dppca';
defaultT         = 50;
defaultTij       = 5;
defaultComm      = 'every';
defaultTHRESH_T  = 1;
defaultTAlpha    = 0.5;
defaultVis       = false;
defaultConvCounter = 4;
defaultConvCounterRMS = 10;

addParameter(p,'InitModel',defaultMODEL);
addParameter(p,'Threshold',defaultTHRESH,@isnumeric);
addParameter(p,'ThresholdA',defaultTHRESHA,@isnumeric);
addParameter(p,'ThresholdR',defaultTHRESHR,@isnumeric);
addParameter(p,'ShowObjPer',defaultITER,@isnumeric);
addParameter(p,'MaxIter',defaultMaxIter);
addParameter(p,'ZeroMean',defaultZeroMean);
addParameter(p,'Eta',defaultETA,@isnumeric);
addParameter(p,'VaryETA',defaultVaryETA);
addParameter(p,'VaryETA_mu',defaultVaryETA_mu,@isnumeric);
addParameter(p,'VaryETA_tau',defaultVaryETA_tau,@isnumeric);
addParameter(p,'VaryETA_max',defaultVaryETA_max,@isnumeric);
addParameter(p,'UseResidual',defaultUseResidual);
addParameter(p,'W_GT',defaultW_GT);
addParameter(p,'FwRes',defaultFwRes);
addParameter(p,'FwResEta',defaultFwResEta,@isnumeric);

addParameter(p,'ModelType',defaultModelType);
addParameter(p,'T',defaultT);
addParameter(p,'Tij',defaultTij);
addParameter(p,'Talpha', defaultTAlpha);
addParameter(p,'Communication', defaultComm);
addParameter(p,'Threshold_T', defaultTHRESH_T);
addParameter(p,'vis', defaultVis);
addParameter(p,'ConvCounter', defaultConvCounter);
addParameter(p,'ConvCounterRMS', defaultConvCounterRMS);

parse(p,varargin{:});

clearvars -global;

% Initialize parameters
model_init = p.Results.InitModel;
THRESH     = p.Results.Threshold;
THRESHA    = p.Results.ThresholdA;
THRESHR    = p.Results.ThresholdR;
iter_obj   = p.Results.ShowObjPer;
COUNTER_MAX = p.Results.MaxIter;
ZeroMean   = p.Results.ZeroMean;
global ETA VaryETA_mu VaryETA_tau VaryETA_max;
ETA        = p.Results.Eta;
VaryETA_mu = p.Results.VaryETA_mu;
VaryETA_tau= p.Results.VaryETA_tau;
VaryETA_max= p.Results.VaryETA_max;
UseResidual = p.Results.UseResidual;
W_GT        = p.Results.W_GT;
FwResEta    = p.Results.FwResEta;
converge_counter_max = p.Results.ConvCounter;
converge_counter_rms_max = p.Results.ConvCounterRMS;

global T Tij;
global THRESH_T;
T          = p.Results.T;
Tij        = p.Results.Tij;
THRESH_T   = p.Results.Threshold_T;

global idxETAUpdateBase idxETAUpdateBase_localObj idxETAUpdateBase_Redusial;
global idxETAUpdateType idxETAUpdateType_homo     idxETAUpdateType_hetero;
global idxETAUpdateStop idxETAUpdateStop_globalT  idxETAUpdateStop_localT;
global idxModelType     idxModelType_AdaptiveETA  idxModelType_localFRS idxModelType_globalFRS;

idxModelType_AdaptiveETA = 1;
idxModelType_localFRS    = 2;
idxModelType_globalFRS   = 3;

idxETAUpdateBase_localObj         = 1;
idxETAUpdateBase_Redusial         = 2;
idxETAUpdateBase_localObjRedusial = 3;

idxETAUpdateType_homo             = 1;
idxETAUpdateType_hetero           = 2;

idxETAUpdateStop_globalT          = 1;
idxETAUpdateStop_localT           = 2;

ModelType  = p.Results.ModelType;
switch ModelType
    case 'dppca_NAP_homo'
        idxModelType     = idxModelType_AdaptiveETA;
        idxETAUpdateBase = idxETAUpdateBase_localObj;
        idxETAUpdateType = idxETAUpdateType_homo;
        idxETAUpdateStop = idxETAUpdateStop_globalT;
    case 'dppca_AP_homo'
        idxModelType     = idxModelType_AdaptiveETA;
        idxETAUpdateBase = idxETAUpdateBase_localObj;
        idxETAUpdateType = idxETAUpdateType_homo;
        idxETAUpdateStop = idxETAUpdateStop_localT;

    case 'dppca_NAP_hetero'
        idxModelType     = idxModelType_AdaptiveETA;
        idxETAUpdateBase = idxETAUpdateBase_localObj;
        idxETAUpdateType = idxETAUpdateType_hetero;
        idxETAUpdateStop = idxETAUpdateStop_globalT;
    case 'dppca_AP_hetero'
        idxModelType     = idxModelType_AdaptiveETA;
        idxETAUpdateBase = idxETAUpdateBase_localObj;
        idxETAUpdateType = idxETAUpdateType_hetero;
        idxETAUpdateStop = idxETAUpdateStop_localT;

    case 'dppca_localVP_homo'
        idxModelType     = idxModelType_AdaptiveETA;
        idxETAUpdateBase = idxETAUpdateBase_Redusial;
        idxETAUpdateType = idxETAUpdateType_homo;
        idxETAUpdateStop = idxETAUpdateStop_globalT;
    case 'dppca_localVPAP_homo'
        idxModelType     = idxModelType_AdaptiveETA;
        idxETAUpdateBase = idxETAUpdateBase_localObjRedusial;
        idxETAUpdateType = idxETAUpdateType_homo;
        idxETAUpdateStop = idxETAUpdateStop_globalT;
    case 'dppca_localVPNAP_homo'
        idxModelType     = idxModelType_AdaptiveETA;
        idxETAUpdateBase = idxETAUpdateBase_localObjRedusial;
        idxETAUpdateType = idxETAUpdateType_homo;
        idxETAUpdateStop = idxETAUpdateStop_localT;

    case 'dppca_localVP_hetero'
        idxModelType     = idxModelType_AdaptiveETA;
        idxETAUpdateBase = idxETAUpdateBase_Redusial;
        idxETAUpdateType = idxETAUpdateType_hetero;
        idxETAUpdateStop = idxETAUpdateStop_globalT;
    case 'dppca_localVPAP_hetero'
        idxModelType     = idxModelType_AdaptiveETA;
        idxETAUpdateBase = idxETAUpdateBase_localObjRedusial;
        idxETAUpdateType = idxETAUpdateType_hetero;
        idxETAUpdateStop = idxETAUpdateStop_globalT;
    case 'dppca_localVPNAP_hetero'
        idxModelType     = idxModelType_AdaptiveETA;
        idxETAUpdateBase = idxETAUpdateBase_localObjRedusial;
        idxETAUpdateType = idxETAUpdateType_hetero;
        idxETAUpdateStop = idxETAUpdateStop_localT;

    case 'dppca_localFRS'
        idxModelType     = idxModelType_localFRS;
    case 'dppca_globalFRS'
        idxModelType     = idxModelType_globalFRS;
    otherwise
        fprintf('Invalid ModelType: %s\n',ModelType);
end

if( idxModelType==idxModelType_AdaptiveETA )
    if( idxETAUpdateStop == idxETAUpdateStop_globalT )
        fn_update_eta_cond = @test_etaupdate_globalT;
    elseif( idxETAUpdateStop == idxETAUpdateStop_localT )
        fn_update_eta_cond = @test_etaupdate_localT;
    end

    if( idxETAUpdateType == idxETAUpdateType_homo )
        if( idxETAUpdateBase == idxETAUpdateBase_localObj )
            fn_update_eta      = @update_eta;
            fn_update_eta_base = @update_eta_initETA_homo;
        elseif( idxETAUpdateBase == idxETAUpdateBase_Redusial )
            fn_update_eta      = @update_eta_Redusial;
            fn_update_eta_base = @update_eta_prevVaryETA_homo;
        elseif( idxETAUpdateBase == idxETAUpdateBase_localObjRedusial )
            fn_update_eta      = @update_eta_Redusial;
            fn_update_eta_base = @update_eta_prevVaryETA_homo;
        end
    elseif( idxETAUpdateType == idxETAUpdateType_hetero )
        if( idxETAUpdateBase == idxETAUpdateBase_localObj )
            fn_update_eta      = @update_eta;
            fn_update_eta_base = @update_eta_initETA_hetero;
        elseif( idxETAUpdateBase == idxETAUpdateBase_Redusial )
            fn_update_eta      = @update_eta_Redusial;
            fn_update_eta_base = @update_eta_prevVaryETA_hetero;
        elseif( idxETAUpdateBase == idxETAUpdateBase_localObjRedusial )
            fn_update_eta      = @update_eta_Redusial;
            fn_update_eta_base = @update_eta_prevVaryETA_hetero;
        end
    end
end

global Talpha;
Talpha      = p.Results.Talpha;
vis        = p.Results.vis;
comm_occasionally = strcmp(p.Results.Communication, 'occasionally');
colors     = hsv(J);

assert(ETA > 0, 'Learning rate (ETA) should be positive!');
assert(VaryETA_mu > 0, 'Please check 0 < VaryETA_mu < 1');
assert(VaryETA_mu < 1, 'Please check 0 < VaryETA_mu < 1');
assert(VaryETA_tau > 0, 'Please check 0 < VaryETA_tau');

% Check validity of initilaization
if (isfield(model_init, 'W') && ~iscell(model_init.W)) || ...
    (isfield(model_init, 'MU') && ~iscell(model_init.MU)) || ...
    (isfield(model_init, 'VAR') && ~iscell(model_init.VAR))
    error('Invalid initialization: please specify distributed model');
end

%--------------------------------------------------------------------------
% We need to broadcast parameters & Lagrange multipliers to neighbors. 
% Here, we use global variables. In real settings, sensors should transmit 
% them over network.
global Wi MUi PRECi oWi oMUi oPRECi;
global LAMBDAi GAMMAi BETAi;

% Local variables
global EZ EZZt;
global Bj MISSj;

% To reduce unnecessary memory copy
global Xj;

%--------------------------------------------------------------------------
% Build Xi for speed up
Xj = cell(J,1);
for i = 1 : J
    Xj{i} = X(:, V == i);
end

% Find i-th node's neighbor set Bi in advance to speed up
Bj = cell(J,1);
for i = 1 : J
    if r == 0
        Bj{i} = [];
    else
        Bj{i} = find(E(i,:) > 0);
    end
end

% Local parameters and auxiliary variables defined here for simplicity.
% In the real environment, these variables reside in each local sensor.
% Initialize parameters
for i = 1 : J
    Wi{i} = model_init.W{i};
    MUi{i} = model_init.MU{i};
    PRECi{i} = 1./model_init.VAR{i};
end

% Initialize Lagrange multipliers. Each edge of each node has a multiplier.
for i = 1 : J
    LAMBDAi{i} = zeros(D, M);
    GAMMAi{i} = zeros(D, 1);
    BETAi{i} = 0;
end

% Initialize latent variables
EZ = cell(J, 1);
EZZt = cell(J, 1);

% Build MISSi for speed up
MISSj = cell(J,1);
for i = 1 : J
    MISSj{i} = isnan(Xj{i});
end

% Adaptive Network Setting
global ETA_ij ETAhalf_ij;
global VaryETA_ij;
global tau_ij sum_tau_ij;
global T_ij T_n_ij;
global Fi Fi_history;

ETA_ij     = zeros(J,J);
ETAhalf_ij = zeros(J,J);
tau_ij     = zeros(J,J);
sum_tau_ij = zeros(J,J);
T_ij       = zeros(J,J);
T_n_ij     = ones(J,J);
for i=1:J
    Bi = Bj{i};
    for j=Bi
        ETA_ij(i,j)     = ETA;
        VaryETA_ij(i,j) = ETA;
        T_ij(i,j)       = Tij;
    end
end
ETAhalf_ij = ETA_ij .* 0.5;

Fi_history     = zeros(COUNTER_MAX,J);
ETA_ij_history = zeros(COUNTER_MAX,J,J);

% residual norms
global pri_norms dual_norms;
pri_norms = cell(J,1);
dual_norms = cell(J,1);

% dual variables
zbar_W = cell(J,1);
zbar_MU = cell(J,1);
zbar_PREC = cell(J,1);
o_zbar_W = cell(J,1);
o_zbar_MU = cell(J,1);
o_zbar_PREC = cell(J,1);

% Learning rate
%ETAhalf = ETA * 0.5;

% Initialize objective function - Lagrangian (we are minimizing this)
oldObjLR = realmax;
objArray = zeros(COUNTER_MAX, J+1); % last one is reserved for total
ssaArray = zeros(COUNTER_MAX,1);
rmsArray = zeros(COUNTER_MAX, 1);
prinormArray = zeros(COUNTER_MAX, 1);
pritolArray = zeros(COUNTER_MAX, 1);
dualnormArray = zeros(COUNTER_MAX, 1);
dualtolArray = zeros(COUNTER_MAX, 1);
if idxModelType == idxModelType_localFRS
    alphaArray = zeros(COUNTER_MAX, J);
    cArray = zeros(COUNTER_MAX, J);
elseif idxModelType == idxModelType_globalFRS
    alphaArray = zeros(COUNTER_MAX, 1);
    cArray = zeros(COUNTER_MAX, 1);
end
Fi = zeros(J, 1);

%--------------------------------------------------------------------------
% Prepare performance measures
converged = 0;
counter = 1;
tic;

% if( vis )
%     figure(1);
%     clf;
% end

if( vis )
    W_history = zeros( [COUNTER_MAX, J, size(Wi{1})] );
end

converge_counter = 0;
converge_counter_rms = 0;

% Main loop
while counter <= COUNTER_MAX

    %----------------------------------------------------------------------
    % Temporarily store parameters to simulate broadcasting and
    % synchronization. All nodes should update before sending their values.
    if ~(exist('didRestart', 'var') && didRestart)
        oWi = Wi;
        oMUi = MUi;
        oPRECi = PRECi;
    end
    
    oLAMBDA = LAMBDAi;
    oGAMMA = GAMMAi;
    oBETA = BETAi;
    
    %----------------------------------------------------------------------
    % In each node: Update parameters locally
    for i = 1 : J
        Fi(i) = dppca_local( M, i, ETA_ij(i,:), ZeroMean );
    end
    
    %----------------------------------------------------------------------
    % In each node: Update Lagrange multipliers
    for i = 1 : J
        Bi = Bj{i};
        for j = Bi
            LAMBDAi{i} = LAMBDAi{i} + ( ETAhalf_ij(i,j) * (Wi{i}    - Wi{j}) );
            GAMMAi{i}  = GAMMAi{i}  + ( ETAhalf_ij(i,j) * (MUi{i}   - MUi{j}) );
            BETAi{i}   = BETAi{i}   + ( ETAhalf_ij(i,j) * (PRECi{i} - PRECi{j}) );
        end
    end
    
    %----------------------------------------------------------------------
    % Stopping criterion checkpoint

    % Compute global objective
    objLR = 0;
    for i = 1 : J
        objLRi = Fi(i);
        Bi = Bj{i};
        for j = Bi
            objLRi = objLRi ...
                    + trace(LAMBDAi{i}' * (Wi{i} - Wi{j})) ...
                    + (GAMMAi{i}' * (MUi{i} - MUi{j})) ...
                    + (BETAi{i} * (PRECi{i} - PRECi{j})) ...
                    + ETAhalf_ij(i,j) * norm(Wi{i} - Wi{j},'fro')^2 ...
                    + ETAhalf_ij(i,j) * norm(MUi{i} - MUi{j},'fro')^2 ...
                    + ETAhalf_ij(i,j) * (PRECi{i} - PRECi{j})^2;
        end
        objArray(counter, i) = objLRi;
        objLR = objLR + objLRi;
    end
    objArray(counter,J+1) = objLR;
    relErr = (objLR - oldObjLR) / abs(oldObjLR);
    oldObjLR = objLR;

    %----------------------------------------------------------------------
    % Compute primal and dual residual (Applies to all parameters)
    % Our constraints: for each edge (i,j), we have W_i=rho_ij, rho_ij=W_j.
    % Primal residual is ||W_i-rho_ij|| + ||rho_ij-W_j||
    % Essentially, this is the sum of constraint violations
    % Dual residual is ETA*(rho_{ij}^{t+1} - rho_{ij}^{t})
    %
    % Compute primal and dual tolerance (Applies to all parameters)
    % Primal tolerance is sqrt(.)*TOL + max{||W_i||, ||rho_{ij}||, ...}*TOL
    % Dual tolerance is sqrt(.)*TOL + ||LAMBDA_i||*TOL ...
    pri_norm = 0;
    dual_norm = 0;
    pri_tol = 0;
    dual_tol = 0;

    if counter == 1
        % initialize dual variables
        for i = 1 : J
            o_zbar_W{i} = zeros(size(Wi{1}));
            o_zbar_MU{i} = zeros(size(MUi{1}));
            o_zbar_PREC{i} = 0;
        end
    else
        % save current dual variables
        o_zbar_W = zbar_W;
        o_zbar_MU = zbar_MU;
        o_zbar_PREC = zbar_PREC;
    end
    
    % compute mean eta value (mean learning rate)
    switch idxModelType
        case idxModelType_AdaptiveETA,
            ETA_i = mean(ETA_ij,2);
        case idxModelType_localFRS,
            ETA_i = repmat(ETA, [1, J]);
        case idxModelType_globalFRS,
            ETA_i = repmat(ETA, [1, J]);
    end

if(0)
    if (idxModelType == idxModelType_localVPAP || idxModelType == idxModelType_localVPNAP)
        for i = 1 : J
            % compute mean and old mean of variables
            zbar_W = Wi{i};
            zbar_MU = MUi{i};
            zbar_PREC = PRECi{i};

            dual_norm_i = 0;
            
            Bi = Bj{i};
            for j = 1 : length(Bi)
                zbar_W = zbar_W + Wi{Bi(j)};
                zbar_MU = zbar_MU + MUi{Bi(j)};
                zbar_PREC = zbar_PREC + PRECi{Bi(j)};

                % compute consensus dual residual
                dual_norm_i = dual_norm_i + norm(ETAhalf_ij(i,Bi(j)) * (Wi{Bi(j)} - oWi{Bi(j)}),'fro');
                dual_norm_i = dual_norm_i + norm(ETAhalf_ij(i,Bi(j)) * (MUi{Bi(j)} - oMUi{Bi(j)}),'fro');
                dual_norm_i = dual_norm_i + norm(ETAhalf_ij(i,Bi(j)) * (PRECi{Bi(j)} - oPRECi{Bi(j)}),'fro');
            end

            dual_norm = dual_norm + dual_norm_i;
            %dual_norm = dual_norm + dual_norm_i / length(Bi);

            zbar_W = zbar_W ./ (length(Bi)+1);
            zbar_MU = zbar_MU ./ (length(Bi)+1);
            zbar_PREC = zbar_PREC ./ (length(Bi)+1);

            % compute consensus primal residual
            pri_norm = pri_norm + norm(Wi{i} - zbar_W,'fro');
            pri_norm = pri_norm + norm(MUi{i} - zbar_MU,'fro');
            pri_norm = pri_norm + norm(PRECi{i} - zbar_PREC,'fro');  
            
            % calculate primal tolerence
            pri_tol_i = 0;
            pri_tol_i = pri_tol_i + THRESHA * sqrt(D*M) + ...
                THRESHR * max(norm(Wi{i},'fro'), norm(zbar_W,'fro'));
            pri_tol_i = pri_tol_i + THRESHA * sqrt(D) + ...
                THRESHR * max(norm(MUi{i},'fro'), norm(zbar_MU,'fro'));
            pri_tol_i = pri_tol_i + THRESHA * sqrt(1) + ...
                THRESHR * max(norm(PRECi{i},'fro'), norm(zbar_PREC,'fro'));            
            pri_tol = pri_tol + pri_tol_i;

            % calculate dual tolerence
            dual_tol_i = 0;
            dual_tol_i = dual_tol_i ...
                + THRESHA * sqrt(D*M) + THRESHR * norm(ETA_i(i) * LAMBDAi{i},'fro');
            dual_tol_i = dual_tol_i ...
                + THRESHA * sqrt(D) + THRESHR * norm(ETA_i(i) * GAMMAi{i},'fro');
            dual_tol_i = dual_tol_i ...
                + THRESHA * sqrt(1) + THRESHR * norm(ETA_i(i) * BETAi{i},'fro');
            dual_tol = dual_tol + dual_tol_i;            
        end
    else
        for i = 1 : J
            % compute mean of neighbors
            zbar_W{i} = zeros(size(Wi{i}));
            zbar_MU{i} = zeros(size(MUi{i}));
            zbar_PREC{i} = 0;

            Bi = Bj{i};
            for j = 1 : length(Bi)
                zbar_W{i} = zbar_W{i} + Wi{Bi(j)};
                zbar_MU{i} = zbar_MU{i} + MUi{Bi(j)};
                zbar_PREC{i} = zbar_PREC{i} + PRECi{Bi(j)};
            end

            zbar_W{i} = zbar_W{i} ./ length(Bi);
            zbar_MU{i} = zbar_MU{i} ./ length(Bi);
            zbar_PREC{i} = zbar_PREC{i} ./ length(Bi);

            % calculate primal tolerence
            pri_tol_i = 0;
            
            param_i = [reshape(Wi{i}, [D*M, 1]); reshape(MUi{i}, [D, 1]); PRECi{i}];
            param_j = [reshape(zbar_W{i}, [D*M, 1]); reshape(zbar_MU{i}, [D, 1]); zbar_PREC{i}];
            THRESHA * sqrt(D*M + D + 1) + THRESHR * max(norm(param_i,'fro'), norm(param_j,'fro'));
            
            
            pri_tol_i = pri_tol_i + THRESHA * sqrt(D*M) + ...
                THRESHR * max(norm(Wi{i},'fro'), norm(zbar_W{i},'fro'));
            pri_tol_i = pri_tol_i + THRESHA * sqrt(D) + ...
                THRESHR * max(norm(MUi{i},'fro'), norm(zbar_MU{i},'fro'));
            pri_tol_i = pri_tol_i + THRESHA * sqrt(1) + ...
                THRESHR * max(norm(PRECi{i},'fro'), norm(zbar_PREC{i},'fro'));            
            pri_tol = pri_tol + pri_tol_i;

            % calculate dual tolerence
            dual_tol_i = 0;
            dual_tol_i = dual_tol_i ...
                + THRESHA * sqrt(D*M) + THRESHR * norm(ETA_i(i) * LAMBDAi{i},'fro');
            dual_tol_i = dual_tol_i ...
                + THRESHA * sqrt(D) + THRESHR * norm(ETA_i(i) * GAMMAi{i},'fro');
            dual_tol_i = dual_tol_i ...
                + THRESHA * sqrt(1) + THRESHR * norm(ETA_i(i) * BETAi{i},'fro');
            dual_tol = dual_tol + dual_tol_i;

            % compute consensus dual residual
            dual_norm_i = 0;
            dual_norm_i = dual_norm_i + ...
                norm(ETA_i(i) * (zbar_W{i} - o_zbar_W{i}),'fro');
            dual_norm_i = dual_norm_i + ...
                norm(ETA_i(i) * (zbar_MU{i} - o_zbar_MU{i}),'fro');
            dual_norm_i = dual_norm_i + ...
                norm(ETA_i(i) * (zbar_PREC{i} - o_zbar_PREC{i}),'fro');
            dual_norm = dual_norm + dual_norm_i;

            % compute consensus primal residual
            pri_norm_i = 0;
            pri_norm_i = pri_norm_i + norm(Wi{i} - zbar_W{i},'fro');
            pri_norm_i = pri_norm_i + norm(MUi{i} - zbar_MU{i},'fro');
            pri_norm_i = pri_norm_i + norm(PRECi{i} - zbar_PREC{i},'fro');
            pri_norm = pri_norm + pri_norm_i;

            % save
            pri_norms{i} = pri_norm_i;
            dual_norms{i} = dual_norm_i;
        end
    end
else
    for i = 1 : J
        % compute mean of neighbors
        zbar_W{i}    = zeros(size(Wi{i}));
        zbar_MU{i}   = zeros(size(MUi{i}));
        zbar_PREC{i} = 0;
        
        Bi = Bj{i};
        for j = 1 : length(Bi)
            zbar_W{i}    = zbar_W{i}    + Wi{Bi(j)};
            zbar_MU{i}   = zbar_MU{i}   + MUi{Bi(j)};
            zbar_PREC{i} = zbar_PREC{i} + PRECi{Bi(j)};
        end

        zbar_W{i}    = zbar_W{i}    ./ length(Bi);
        zbar_MU{i}   = zbar_MU{i}   ./ length(Bi);
        zbar_PREC{i} = zbar_PREC{i} ./ length(Bi);

        % calculate primal tolerence
        pri_tol_norm = max( sqrt(   norm(Wi{i},       'fro')^2          ...
                                  + norm(MUi{i},      'fro')^2          ...
                                  + norm(PRECi{i},    'fro')^2     ),   ...
                            sqrt(   norm(zbar_W{i},   'fro')^2          ...
                                  + norm(zbar_MU{i},  'fro')^2          ...
                                  + norm(zbar_PREC{i},'fro')^2     ) );
        pri_tol_i = THRESHA*sqrt(D*M + D + 1) + THRESHR*pri_tol_norm;
        pri_tol = pri_tol + pri_tol_i;

        % calculate dual tolerence
        dual_tol_norm = sqrt(   norm(ETA_i(i) * LAMBDAi{i},'fro')^2 ...
                              + norm(ETA_i(i) * GAMMAi{i}, 'fro')^2 ...
                              + norm(ETA_i(i) * BETAi{i},  'fro')^2     );
        dual_tol_i = THRESHA*sqrt(D*M + D + 1) + THRESHR*dual_tol_norm;
        dual_tol = dual_tol + dual_tol_i;

        % compute consensus dual residual
        dual_norm_i = sqrt(  norm(ETA_i(i) * (zbar_W{i}    - o_zbar_W{i}),   'fro')^2 ...
                           + norm(ETA_i(i) * (zbar_MU{i}   - o_zbar_MU{i}),  'fro')^2 ...
                           + norm(ETA_i(i) * (zbar_PREC{i} - o_zbar_PREC{i}),'fro')^2    );
        dual_norm = dual_norm + dual_norm_i;
        
        % compute consensus primal residual
        pri_norm_i = sqrt(   norm(Wi{i}    - zbar_W{i},   'fro')^2 ...
                           + norm(MUi{i}   - zbar_MU{i},  'fro')^2 ...
                           + norm(PRECi{i} - zbar_PREC{i},'fro')^2     );
        pri_norm = pri_norm + pri_norm_i;

        % save
        pri_norms{i}  = pri_norm_i;
        dual_norms{i} = dual_norm_i;
    end

    prinormArray(counter)  = pri_norm;
    pritolArray(counter)   = pri_tol;
    dualnormArray(counter) = dual_norm;
    dualtolArray(counter)  = dual_tol;

end
    % Compute c_k (centralized)
    if (idxModelType == idxModelType_localFRS)
        % compute c_k
        show_iter = false;
        % for each node...
        for i = 1 : J
            term1 = 0;
            term1 = term1 + norm(LAMBDAi{i} - oLAMBDA{i}, 'fro');
            term1 = term1 + norm(GAMMAi{i} - oGAMMA{i}, 'fro');
            term1 = term1 + norm(BETAi{i} - oBETA{i}, 'fro');
            term2 = 0;
            term2 = term2 + norm(zbar_W{i} - o_zbar_W{i}, 'fro');
            term2 = term2 + norm(zbar_MU{i} - o_zbar_MU{i}, 'fro');
            term2 = term2 + norm(zbar_PREC{i} - o_zbar_PREC{i}, 'fro');
            cArray(counter,i) = (2/ETA) * term1 + (2*ETA) * term2;
            if counter > 1
                if cArray(counter,i) < cArray(counter-1,i) * FwResEta
                    alphaArray(counter,i) = ( 1+sqrt(1+4*alphaArray(counter-1,i).^2) )/2;
                    alphaWeight = (alphaArray(counter-1,i) - 1) / alphaArray(counter,i);
                    zbar_W{i} = zbar_W{i}+ alphaWeight * (zbar_W{i} - o_zbar_W{i});
                    zbar_MU{i} = zbar_MU{i}+ alphaWeight * (zbar_MU{i} - o_zbar_MU{i});
                    zbar_PREC{i} = zbar_PREC{i} + alphaWeight * (zbar_PREC{i} - o_zbar_PREC{i});
                    LAMBDAi{i} = LAMBDAi{i} + alphaWeight * (LAMBDAi{i} - oLAMBDA{i});
                    GAMMAi{i} = GAMMAi{i} + alphaWeight * (GAMMAi{i} - oGAMMA{i});
                    BETAi{i} = BETAi{i} + alphaWeight * (BETAi{i} - oBETA{i});
                    didRestart = false;
                else
                    % We can't ensure monotonicity here...simply ignore the
                    % acceleration if the decrease is smaller than FwResEta
                    alphaArray(counter,i) = 1;
                    cArray(counter,i) = cArray(counter-1,i) / FwResEta;
                    zbar_W{i} = o_zbar_W{i};
                    zbar_MU{i} = o_zbar_MU{i};
                    zbar_PREC{i} = o_zbar_PREC{i};
                    LAMBDAi{i} = oLAMBDA{i};
                    GAMMAi{i} = oGAMMA{i};
                    BETAi{i} = oBETA{i};
                    didRestart = true;
                end
                if iter_obj > 0 && mod(counter, iter_obj) == 0 && (show_iter > 1)
                    if didRestart
                        fprintf('  Node %d, c_{k} = %.2f, c_{k-1} = %.2f, alpha_{k} = %.2f, alpha_{k-1} = %.2f (RESTARTED)\n', ...
                            i, cArray(counter,i), cArray(counter-1,i), alphaArray(counter,i), alphaArray(counter-1,i));
                    else
                        fprintf('  Node %d, c_{k} = %.2f, c_{k-1} = %.2f, alpha_{k} = %.2f, alpha_{k-1} = %.2f\n', ...
                            i, cArray(counter,i), cArray(counter-1,i), alphaArray(counter,i), alphaArray(counter-1,i));
                    end
                end
            else
                alphaArray(counter,i) = (1 + sqrt(5)) / 2;
                if show_iter > 1
                    fprintf('  Node %d, c_{k} = %.2f, c_{k-1} = n/a, alpha_{k} = %.2f, alpha_{k-1} = n/a (INIT)\n', ...
                        i, cArray(counter,i), alphaArray(counter,i));
                end
            end
            if show_iter && counter > 1
                fprintf('  Node %d, c_{k} = %.2f, c_{k-1} = %.2f, relerr = %f\n', ...
                    i, cArray(counter,i), cArray(counter-1,i), ...
                    abs(cArray(counter,i) - cArray(counter-1,i)) / cArray(counter-1,i));
            end
        end
    elseif (idxModelType == idxModelType_globalFRS)
        % compute c_k
        show_iter = false;
        term1 = 0;
        term2 = 0;
        for i = 1 : J
            term1 = term1 + norm(LAMBDAi{i} - oLAMBDA{i}, 'fro');
            term1 = term1 + norm(GAMMAi{i} - oGAMMA{i}, 'fro');
            term1 = term1 + norm(BETAi{i} - oBETA{i}, 'fro');
            term2 = term2 + norm(zbar_W{i} - o_zbar_W{i}, 'fro');
            term2 = term2 + norm(zbar_MU{i} - o_zbar_MU{i}, 'fro');
            term2 = term2 + norm(zbar_PREC{i} - o_zbar_PREC{i}, 'fro');
        end
        cArray(counter) = (2/ETA) * term1 + (2*ETA) * term2;
        if counter > 1
            if cArray(counter) < cArray(counter-1) * FwResEta
                alphaArray(counter) = ( 1+sqrt(1+4*alphaArray(counter-1).^2) )/2;
                alphaWeight = (alphaArray(counter-1) - 1) / alphaArray(counter);
                for i = 1 : J
                    zbar_W{i} = zbar_W{i}+ alphaWeight * (zbar_W{i} - o_zbar_W{i});
                    zbar_MU{i} = zbar_MU{i}+ alphaWeight * (zbar_MU{i} - o_zbar_MU{i});
                    zbar_PREC{i} = zbar_PREC{i} + alphaWeight * (zbar_PREC{i} - o_zbar_PREC{i});
                    LAMBDAi{i} = LAMBDAi{i} + alphaWeight * (LAMBDAi{i} - oLAMBDA{i});
                    GAMMAi{i} = GAMMAi{i} + alphaWeight * (GAMMAi{i} - oGAMMA{i});
                    BETAi{i} = BETAi{i} + alphaWeight * (BETAi{i} - oBETA{i});
                end
                didRestart = false;
            else
                % We can't ensure monotonicity here...simply ignore the
                % acceleration if the decrease is smaller than FwResEta
                alphaArray(counter) = 1;
                cArray(counter) = cArray(counter-1) / FwResEta;
                for i = 1 : J
                    zbar_W{i} = o_zbar_W{i};
                    zbar_MU{i} = o_zbar_MU{i};
                    zbar_PREC{i} = o_zbar_PREC{i};
                    LAMBDAi{i} = oLAMBDA{i};
                    GAMMAi{i} = oGAMMA{i};
                    BETAi{i} = oBETA{i};
                end
                didRestart = true;
            end
            if iter_obj > 0 && mod(counter, iter_obj) == 0 && (show_iter > 1)
                if didRestart
                    fprintf('  c_{k} = %.2f, c_{k-1} = %.2f, alpha_{k} = %.2f, alpha_{k-1} = %.2f (RESTARTED)\n', ...
                        cArray(counter), cArray(counter-1), alphaArray(counter), alphaArray(counter-1));
                else
                    fprintf('  c_{k} = %.2f, c_{k-1} = %.2f, alpha_{k} = %.2f, alpha_{k-1} = %.2f\n', ...
                        cArray(counter), cArray(counter-1), alphaArray(counter), alphaArray(counter-1));
                end
            end
        else
            alphaArray(counter) = (1 + sqrt(5)) / 2;
            if show_iter > 1
                fprintf('  c_{k} = %.2f, c_{k-1} = n/a, alpha_{k} = %.2f, alpha_{k-1} = n/a (INIT)\n', ...
                    cArray(counter), alphaArray(counter));
            end
        end
    end
    
    % Calculate subspace angle if requested
    if ~isempty(W_GT)
        tmpW = Wi;
        for i = 1 : J
            tmpW{i} = remove_rot_amb_temp(Wi{i}, EZ{i}, EZZt{i}, Xj{i});
        end
        ssaArray(counter) = calc_ppca_max_ssa_gt(tmpW, W_GT, J);
    end    
    % Calculate reprojection error
    rmsArray(counter) = calc_ppca_rms(X, Wi, EZ, MUi);
    
    % Show progress if requested
    if iter_obj > 0 && mod(counter, iter_obj) == 0
        if ~UseResidual
            fprintf('Iter = %d:  Cost = %f (rel %3.2f%%), RMS = %.2e (J = %d, ETA = %f)\n', ...
                counter, objLR, relErr*100, ...
                calc_ppca_rms(Xj, Wi, EZ, MUi), ...
                J, ETA);
%                 calc_ppca_max_ssa(oWi, Wi), ...
        else
            fprintf('Iter = %d:  p_norm = %f / p_tol = %f / d_norm = %f / d_tol = %f, RMS = %.2e (J = %d, ETA = %f)\n', ...
                counter, pri_norm, pri_tol, dual_norm, dual_tol, ...
                calc_ppca_rms(Xj, Wi, EZ, MUi), ...
                J, ETA);
%                 calc_ppca_max_ssa(oWi, Wi), ...
        end
    end

%     if iter_obj > 0 && mod(counter, iter_obj) == 0
%         fprintf('--Iter = %d:  Cost = %f (rel %3.3f%%)\n',counter, objLR, relErr*100);
%     end
    
    % Check whether it has converged
    if (~UseResidual && idxModelType == idxModelType_localFRS)
        if abs(relErr) < THRESH
            converge_counter = converge_counter + 1;
            if converge_counter >= converge_counter_max
                converged = 1;
                break;
            end
        else
            converged = 1;
            for i = 1 : J
                % compute relative change of c_k
                newC = cArray(counter,i);
                if counter <= 1
                    oldC = realmax;
                else
                    oldC = cArray(counter-1,i);
                end
                if (abs(newC - oldC)/oldC) > THRESH
                    converged = 0;
                    break;
                end
            end
            if converged
                break;
            end
        end
    elseif (~UseResidual && idxModelType == idxModelType_globalFRS)
        newC = cArray(counter);
        if counter <= 1
            oldC = realmax;
        else
            oldC = cArray(counter-1);
        end
        if abs(relErr) < THRESH && (abs(newC - oldC) / oldC) < THRESH
            converge_counter = converge_counter + 1;
            if converge_counter >= converge_counter_max
                converged = 1;
                break;
            end
        end
    elseif UseResidual
%         if J == 1 && abs(relErr) < THRESH
%             converge_counter = converge_counter + 1;                        
        if J > 1 && pri_norm <= pri_tol && dual_norm <= dual_tol
%             converged = 1;
%             break;
%         else
%             converged = 0;
            converge_counter = converge_counter + 1;                        
        else
            converge_counter = 0;
        end
        
        if converge_counter >= converge_counter_max
            converged = 1;
            break;
        end
    else
        if abs(relErr) < THRESH
            converge_counter = converge_counter + 1;
            if converge_counter >= converge_counter_max
                converged = 1;
                break;
            end
        else
            converge_counter = 0;    
        end
    end
    
    % if reprojection error increased after washing out initial effects,
    % stop the iteration. In the worst case, it takes (J-1) iterations for
    % all nodes in the network to receive information of all the other
    % nodes at least once (chain). Thus we take J*2 as the threshold.
    if J < 10
        rmsbound = J * 2;
    else
        rmsbound = J;
    end    
    if counter > rmsbound && rmsArray(counter-1) < rmsArray(counter)
        converge_counter_rms = converge_counter_rms + 1;
        if converge_counter_rms >= converge_counter_rms_max
            converged = 1;
            break;
%         elseif converge_counter_rms > 0
%             converge_counter_rms = converge_counter_rms - 1;
        end
    end
    
    if( vis )
        for i=1:J
            W_history(counter,i,:,:) = Wi{i};
        end
    end
    
    %----------------------------------------------------------------------
    % In each node: Update tau_ij
    if( (idxETAUpdateBase == idxETAUpdateBase_localObj        ) || ...
        (idxETAUpdateBase == idxETAUpdateBase_localObjRedusial)   )
        for i = 1 : J
            update_tau(i, counter);
        end
    end

    %----------------------------------------------------------------------
    % In each node: Update T_ij
    if( (idxETAUpdateStop == idxETAUpdateStop_localT) )
        for i = 1 : J
            update_T(i, counter);
        end
    end
    
    %----------------------------------------------------------------------
    % In each node: Update ETA_ij
    for i = 1:J
        fn_update_eta(i,counter,fn_update_eta_cond, fn_update_eta_base);
    end
    ETAhalf_ij = ETA_ij .* 0.5;

    if(vis)
        hold on;
        clf;

        for i = 1:J
            subplot(J,1,i);
            hold on;

            for j=Bj{i}
                ETA_ij_history(counter,i,j) = ETA_ij(i,j);                
                plot(1:counter,ETA_ij_history(1:counter,i,j),'color',colors(j,:));
            end        
        end   

        drawnow;
    end
        
    % Increase counter
    counter = counter + 1;
end

% Check convergence
if converged ~= 1
    counter = COUNTER_MAX;
    %fprintf('Could not converge within %d iterations.\n', COUNTER_MAX);
end

% Compute performance measures
eTIME = toc;
eITER = counter;

% Fill in remaining objective function value slots.
if counter < COUNTER_MAX
    % fill in the remaining objArray
    for i = 1:J+1
        objArray(counter+1:COUNTER_MAX,i) = ...
            repmat(objArray(counter,i), [COUNTER_MAX - counter, 1]);
    end
end

% Truncate to save memory
ssaArray = ssaArray(1:counter);
rmsArray = rmsArray(1:counter);
prinormArray  = prinormArray(1:counter);
pritolArray   = pritolArray(1:counter);
dualnormArray = dualnormArray(1:counter);
dualtolArray  = dualtolArray(1:counter);


% Remove rotational ambiguity
% for idx = 1 : J
%     [~, ~, V] = svd(Wi{idx});
%     Wi{idx} = Wi{idx} * V;
%     EZ{idx} = V'*EZ{idx};
%     for idy = 1 : size(Xj{idx},2)
%         EZZt{idx}(:,:,idy) = V' * EZZt{idx}(:,:,idy) * V;
%     end
% 
%     % OR...
%     [t1, t2] = remove_rot_amb(EZ{idx}', Wi(:,:,idx)');
%     EZ{idx} = t1';
%     Wi(:,:,idx) = t2';%
% 
%     % Correct rotational ambiguity (Sec. 4.1, 4.2, Ilin and Raiko, 2010)
%     [D, Ni] = size(Xj{idx});
% 
%     mEZn = mean(EZ{idx},2);
%     MUi(:,idx) = MUi(:,idx) + Wi(:,:,idx) * mEZn;
%     EZ{idx} = bsxfun(@minus, EZ{idx}, mEZn);
%     
%     vEZnZnt = EZ{idx} * EZ{idx}';
%     for idn = 1:Ni
%         vEZnZnt = vEZnZnt + EZZt{idx}(:,:,idn);
%     end
%     vEZnZnt = vEZnZnt ./ Ni;
%     [Uz, Dz] = eig(vEZnZnt);
%     Dz = sqrt(Dz);
%     vWWt = ((Wi(:,:,idx) * Uz * Dz)' * (Wi(:,:,idx) * Uz * Dz)) ./ D;
%     [Vw, Dw] = eig(vWWt);
%     [~, I] = sort(-diag(Dw));
%     Vw = Vw(:,I);
% 
%     Wi(:,:,idx) = Wi(:,:,idx) * Uz * Dz * Vw;
%     R = Vw' * diag(1 ./ diag(Dz)) * Uz';
%     EZ{idx} = R * EZ{idx};
%     for idn = 1:Ni
%         EZZt{idx}(:,:,idn) = R * EZZt{idx}(:,:,idn) * R';
%     end
% end

% Assign return values
W = cell(J,1);
MU = cell(J,1);
VAR = cell(J,1);
for i = 1 : J
    W{i} = Wi{i};
    MU{i} = MUi{i};
    VAR{i} = 1./PRECi{i};
end

ETA_ij_history = ETA_ij_history(1:counter,:,:);
%tau_ij_history = tau_ij_history(1:counter,:,:);
%T_ij_history   = T_ij_history(1:counter,:,:);
%comm_histroy   = comm_histroy(1:counter,:,:);

% Create structure
if( vis )
model = structure( ...
    W, MU, VAR, ...
    EZ, EZZt, ...
    eITER, eTIME, objArray, ssaArray, rmsArray, ...
    LAMBDAi, GAMMAi, BETAi, ...
    ETA_ij_history, W_history, ...
    prinormArray, pritolArray, dualnormArray, dualtolArray );
    %ETA_ij_history, tau_ij_history, T_ij_history, comm_histroy, W_history);

else
    
model = structure( ...
    W, MU, VAR, ...
    EZ, EZZt, ...
    eITER, eTIME, objArray, ssaArray, rmsArray, ...
    LAMBDAi, GAMMAi, BETAi, ...
    prinormArray, pritolArray, dualnormArray, dualtolArray );

end

% Clean up
clearvars -global Xj;
clearvars -global Wi MUi PRECi oWi oMUi oPRECi;
clearvars -global EZ EZZt;
clearvars -global Bj MISSj;
clearvars -global LAMBDAi GAMMAi BETAi;

clearvars -except model;

end

function [tmpW] = remove_rot_amb_temp(tW, tEZ, tEZZt, tXj)

[D, Ni] = size(tXj);

mEZn = mean(tEZ,2);
tEZ = bsxfun(@minus, tEZ, mEZn);

vEZnZnt = tEZ * tEZ';
for idn = 1:Ni
    vEZnZnt = vEZnZnt + tEZZt(:,:,idn);
end
vEZnZnt = vEZnZnt ./ Ni;
[Uz, Dz] = eig(vEZnZnt);
Dz = sqrt(Dz);
vWWt = ((tW * Uz * Dz)' * (tW * Uz * Dz)) ./ D;
[Vw, Dw] = eig(vWWt);
[~, I] = sort(-diag(Dw));
Vw = Vw(:,I);

tmpW = tW * Uz * Dz * Vw;

end

function [F_new] = dppca_local( M, i, ETAj, isMeanZero )
% DPPCA_LOCAL  D-PPCA Local Update
% 
% Input
%  M     : Projected dimension
%  i     : Current node index
%  ETA   : Scalar Learning ratio
%  isMeanZero : True if we don't update MUi
%
% Output
%  F_new    : 1 x 1 scalar computed optimization forumla (first term only)

% Parameters and latent space variables: Not that only the three model 
% parameters need to be transmitted. Other variables are defined as global
% just for simple and easy-to-understand implementation.
global Wi MUi PRECi oWi oMUi oPRECi;
global LAMBDAi GAMMAi BETAi;
global EZ EZZt;
global Bj MISSj;
global Xj;

% Take i-th node
Xi = Xj{i};
Bi = Bj{i};
MISSi = MISSj{i};

% Get size of this samples and ball of this node
[D, Ni] = size(Xi);
cBj = length(Bi);

% Sum of ETAj
sumETAj = sum(ETAj(Bi));

% Initialize latent variables (for loop implementation)
EZn = zeros(M, Ni);
EZnZnt = zeros(M, M, Ni);

%--------------------------------------------------------------------------
% E-step

for n = 1 : Ni
    % Get indicies of available features
    DcI = (MISSi(:,n) == 0);    
    
    % Compute Mi = Wi'Wi + VARi*I first
    Wc = Wi{i}(DcI,:);
    Mi = Wc' * Wc + 1/PRECi{i} * eye(M);

    % E[Zn] = Mi^(-1) * Wi' * (Xin - MUi)
    % Currently M x N
    EZn(:,n) = Mi \ Wc' * (Xi(DcI,n) - MUi{i}(DcI));

    % E[z_n z_n'] = VAR * Mi^(-1) + E[z_n]E[z_n]'
    % Currently M x M
    EZnZnt(:,:,n) = inv(Mi) * (1/PRECi{i}) + EZn(:,n) * EZn(:,n)';
end

%--------------------------------------------------------------------------
% M-step

% Update Wi
W_new = zeros(D, M);
% One dimension at a time
for d = 1 : D
    % Get non-missing point indexes
    NcI = (MISSi(d,:) == 0);

    %W_new1 = sum( EZnZnt(:,:,NcI), 3 ) * oPRECi{i} + 2*ETA*cBj*eye(M);
    W_new1 = sum( EZnZnt(:,:,NcI), 3 ) * oPRECi{i} + 2*sumETAj*eye(M);
    W_new2 = ( (Xi(d,NcI) - oMUi{i}(d)) * EZn(:,NcI)' ) * oPRECi{i};
    W_new3 = 2 * LAMBDAi{i}(d,:);
    W_new4 = zeros(1, M);
    
%     for j = 1:cBj
%         W_new4 = W_new4 + (oWi{i}(d,:) + oWi{Bi(j)}(d,:));
%     end
%     
%     % Update
%     if sum(sum(W_new1)) < eps
%         W_new(d,:) = zeros(1,M);
%     else
%         W_new(d,:) = (W_new2 - W_new3 + ETA * W_new4) / W_new1;
%     end
    for j = Bi
        W_new4 = W_new4 + ETAj(j)*(oWi{i}(d,:) + oWi{j}(d,:));
    end
    
    % Update
    if sum(sum(W_new1)) < eps
        W_new(d,:) = zeros(1,M);
    else
        W_new(d,:) = (W_new2 - W_new3 + W_new4) / W_new1;
    end
end

% Update MUi
MU_new = zeros(D,1);
if ~isMeanZero
    for d = 1 : D
        % Get non-missing point indexes
        NcI = (MISSi(d,:) == 0);
        
        %MU_new1 = sum(NcI) * oPRECi{i} + 2*ETA*cBj;
        MU_new1 = sum(NcI) * oPRECi{i} + 2*sumETAj;
        MU_new2 = oPRECi{i} * sum( Xi(d,NcI) - W_new(d,:) * EZn(:,NcI), 2 );
        MU_new3 = 2 * GAMMAi{i}(d);
        MU_new4 = 0;
%         for j = 1 : cBj
%             MU_new4 = MU_new4 + ( oMUi{i}(d) + oMUi{Bi(j)}(d) );
%         end        
%         % Update
%         MU_new(d) = (MU_new2 - MU_new3 + ETA * MU_new4) / MU_new1;
        for j = Bi
            MU_new4 = MU_new4 + (ETAj(j))*( oMUi{i}(d) + oMUi{j}(d) );
        end        
        % Update
        MU_new(d) = (MU_new2 - MU_new3 + MU_new4) / MU_new1;
    end
end

% Update PRECi (by solve for VARi^(-1))
%PREC_new1 = 2 * ETA * cBj;
PREC_new1 = 2 * sumETAj;
PREC_new21 = 2 * BETAi{i};
PREC_new22 = 0;
for j = Bi
    %PREC_new22 = PREC_new22 + ETA * (oPRECi{i} + oPRECi{Bi(j)});
    PREC_new22 = PREC_new22 + ETAj(j) * (oPRECi{i} + oPRECi{j});
end
PREC_new23 = 0;
PREC_new24 = 0;
PREC_new4 = 0;
for n = 1:Ni
    % Get indices of available features
    DcI = (MISSi(:,n) == 0);
    PREC_new4 = PREC_new4 + sum(DcI);
    Wc = W_new(DcI,:);
    
    PREC_new23 = PREC_new23 + EZn(:,n)' * Wc' * (Xi(DcI,n) - MU_new(DcI));
    PREC_new24 = PREC_new24 + 0.5 * ( norm( Xi(DcI,n) - MU_new(DcI), 2 )^2 ...
        + trace( EZnZnt(:,:,n) * (Wc' * Wc) ) );
end
PREC_new2 = PREC_new21 - PREC_new22 - PREC_new23 + PREC_new24;
PREC_new3 = -PREC_new4 / 2;
PREC_new = roots([PREC_new1, PREC_new2, PREC_new3]);

% We follow larger, real solution.
if length(PREC_new) > 1
    PREC_new = max(PREC_new);
end
if abs(imag(PREC_new)) ~= 0i
    error('No real solution!');
    % We shouldn't reach here since both solutions are not real...
end
if PREC_new < 0
    error('Negative precicion!');
end

% Correct rotational ambiguity (Sec. 4.1, 4.2, Ilin and Raiko, 2010)
% mEZn = mean(EZn,2);
% MU_new = MU_new + W_new * mEZn;
% EZn = bsxfun(@minus, EZn, mEZn);
% 
% vEZnZnt = EZn * EZn';
% for idn = 1:Ni
%     vEZnZnt = vEZnZnt + EZnZnt(:,:,idn);
% end
% vEZnZnt = vEZnZnt ./ Ni;
% [Uz, Dz] = eig(vEZnZnt);
% Dz = sqrt(Dz);
% vWWt = ((W_new * Uz * Dz)' * (W_new * Uz * Dz)) ./ D;
% [Vw, Dw] = eig(vWWt);
% [~, I] = sort(-diag(Dw));
% Vw = Vw(:,I);
% 
% W_new = W_new * Uz * Dz * Vw;
% R = Vw' * diag(1 ./ diag(Dz)) * Uz';
% EZn = R * EZn;
% for idn = 1:Ni
%     EZnZnt(:,:,idn) = R * EZnZnt(:,:,idn) * R';
% end

% Update parameter values
Wi{i} = W_new;
MUi{i} = MU_new;
PRECi{i} = PREC_new;

EZ{i} = EZn;
EZZt{i} = EZnZnt;

F_new = f_i( i, W_new, MU_new, PREC_new );

end

function F_new = f_i(i, W, MU, PREC)
global Xj;
global MISSj;
global EZ;
global EZZt;

[~, Ni] = size(Xj{i});

% Compute data log likelihood (we don't need to compute constant)
obj_val1 = 0;
obj_val2 = 0;
obj_val3 = 0;
obj_val4 = 0;
obj_val5 = 0;
for n = 1:Ni
    Id = (MISSj{i}(:,n) == 0);
    Xc = Xj{i}(Id,n);
    MUc = MU(Id);
    Wc = W(Id,:);        

    obj_val1 = obj_val1 + 0.5 * trace(EZZt{i}(:,:,n));
    obj_val2 = obj_val2 + 0.5 * sum(Id) * log(2 * pi * PREC);
    obj_val3 = obj_val3 + (PREC/2) * norm(Xc - MUc, 2).^2;
    obj_val4 = obj_val4 + (PREC/2) * trace(EZZt{i}(:,:,n) * (Wc' * Wc));
    obj_val5 = obj_val5 + PREC * EZ{i}(:,n)' * Wc' * (Xc - MUc);
end
F_new = obj_val1 - obj_val2 + obj_val3 + obj_val4 - obj_val5;

end

function update_tau(i, counter)
global Fi Bj;
global tau_ij sum_tau_ij;
global T_ij;
global Wi MUi PRECi;

%if( counter==1 )
%    return;
%end 

b_continue = false;
for j=Bj{i}
if( sum_tau_ij(i,j) < T_ij(i,j) )
    b_continue = true;
end
end
if( b_continue==false )
return;
end

cBi = length(Bj{i});

Fi_j = zeros(cBi,1);
Fi_i = Fi(i);

f_i_max = Fi_i;
f_i_min = Fi_i;
for idx_Bi=1:cBi
    j = Bj{i}(idx_Bi);
    f = f_i( i, Wi{j}, MUi{j}, PRECi{j} );

    Fi_j(idx_Bi) = f;
    if( f_i_max < f )
        f_i_max = f;
    end
    if( f_i_min > f )
        f_i_min = f;
    end
end

fn_i = (Fi_i - f_i_min)/(f_i_max - f_i_min) + 1;
fn_j = zeros(cBi,1);
for idx_Bi=1:cBi
    j = Bj{i}(idx_Bi);

    %if( f_i_max == Fi_i )
    %    tau_ij(i,j) = 0;
    %else
        fn_j(idx_Bi) = (Fi_j(idx_Bi) - f_i_min)/(f_i_max - f_i_min) + 1;
        tau_ij(i,j) = (fn_i/fn_j(idx_Bi)) - 1;
        
        %tau_ij(i,j) = tau_ij(i,j) * tau_ij(i,j);
    %end

    sum_tau_ij(i,j) = sum_tau_ij(i,j) + abs(tau_ij(i,j));
end

end

function update_T(i, counter)
global Bj;
global Fi_history;
global sum_tau_ij;
global T_ij T_n_ij Tij THRESH_T;
global Talpha;

if( counter < 2 )
    return;
end

Fi_i_cur  = Fi_history(counter,i);
Fi_i_prev = Fi_history(counter-1,i);

if( abs((Fi_i_cur - Fi_i_prev)) > THRESH_T )
    for j=Bj{i}
        if( sum_tau_ij(i,j) >= T_ij(i,j) )
            T_ij(i,j) = T_ij(i,j) + Tij * (Talpha^(T_n_ij(i,j)));
            T_n_ij(i,j) = T_n_ij(i,j) + 1;
        end
    end
end
end

function update_eta_Redusial(i, counter, fn_update_eta_cond, fn_update_eta_base)
global Bj;
global ETA_ij;
global tau_ij sum_tau_ij;
global T_ij;
global VaryETA_mu VaryETA_tau VaryETA_ij;
global pri_norms dual_norms;

pri_norm = pri_norms{i};
dual_norm = dual_norms{i};

if pri_norm * VaryETA_mu >= dual_norm
    for j=Bj{i}
        if( fn_update_eta_cond(i,j,counter) )
            VaryETA_ij(i,j) = VaryETA_ij(i,j) * (1 + VaryETA_tau);
            sum_tau_ij(i,j) = sum_tau_ij(i,j) + VaryETA_tau;
        end
    end
elseif dual_norm * VaryETA_mu >= pri_norm
    for j=Bj{i}
        if( fn_update_eta_cond(i,j,counter) )
            VaryETA_ij(i,j) = VaryETA_ij(i,j) / (1 + VaryETA_tau);
            sum_tau_ij(i,j) = sum_tau_ij(i,j) + VaryETA_tau;
        end
    end
end

update_eta(i,counter, fn_update_eta_cond, fn_update_eta_base);

end

function update_eta(i, counter, fn_update_eta_cond, fn_update_eta_base)
global Bj;
for j=Bj{i}
    fn_update_eta_base(i,j,counter,fn_update_eta_cond);
end
end

function update_eta_prevVaryETA_homo(i,j,counter,fn_test_cond)
global VaryETA_ij;
global ETA ETA_ij;
global tau_ij;
if( fn_test_cond(i,j,counter) )
    ETA_ij(i,j) = VaryETA_ij(i,j) * (1 + tau_ij(i,j));	
else
    ETA_ij(i,j) = ETA;
end
end

function update_eta_prevVaryETA_hetero(i,j,counter,fn_test_cond)
global VaryETA_ij;
global ETA ETA_ij;
global tau_ij;
if( fn_test_cond(i,j,counter) )
    ETA_ij(i,j) = VaryETA_ij(i,j) * (1 + tau_ij(i,j));	
else
    ETA_ij(i,j) = VaryETA_ij(i,j);    
end
end

function update_eta_prevETA_homo(i,j,counter,fn_test_cond)
global ETA ETA_ij;
global tau_ij;
if( fn_test_cond(i,j,counter) )
    ETA_ij(i,j) = ETA_ij(i,j) * (1 + tau_ij(i,j));	
else
    ETA_ij(i,j) = ETA;    
end
end

function update_eta_prevETA_hetero(i,j,counter,fn_test_cond)
global ETA ETA_ij;
global tau_ij;
if( fn_test_cond(i,j,counter) )
    ETA_ij(i,j) = ETA_ij(i,j) * (1 + tau_ij(i,j));	
else
    ETA_ij(i,j) = ETA_ij(i,j);    
end
end

function update_eta_initETA_hetero(i,j,counter,fn_test_cond)
global ETA ETA_ij;
global tau_ij;
if( fn_test_cond(i,j,counter) )
    ETA_ij(i,j) = ETA * (1 + tau_ij(i,j));	
else
    ETA_ij(i,j) = ETA_ij(i,j);
end
end

function update_eta_initETA_homo(i,j,counter,fn_test_cond)
global ETA ETA_ij;
global tau_ij;
if( fn_test_cond(i,j,counter) )
    ETA_ij(i,j) = ETA * (1 + tau_ij(i,j));	
else
    ETA_ij(i,j) = ETA;    
end
end

function do_update = test_etaupdate_globalT(i,j,counter)
global T;
    do_update = ( counter < T );
end

function do_update = test_etaupdate_localT(i,j,counter)
global sum_tau_ij T_ij;
    do_update = ( sum_tau_ij(i,j) < T_ij(i,j) );
end
