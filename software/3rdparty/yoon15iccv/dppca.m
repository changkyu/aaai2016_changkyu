function model = dppca( X, M, V, E, varargin )
% DPPCA      Distributed Probablistic PCA (D-PPCA)
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
%  MaxIter    : Maximum iterations (Def: 1000)
%  ZeroMean   : True if we enforce the mean to be zero. (Def: false)
%  Eta        : Scalar of learning rate (Def: 10)
%  VaryETA    : Varying Eta value as in (He, et al. 2000)
%  VaryETA_mu : Parameter mu  for (He, et al. 2000) (Def: 0.1) (0 < mu < 1)
%  VaryETA_tau: Parameter tau for (He, et al. 2000) (Def: 1)   (0 < tau)
%  VaryETA_max: Max iterations    (He, et al. 2000) (Def: 50)
%  UseResidual: Use residual for stopping criterion (Def: false)
%  W_GT       : Ground truth projection matrix if you want iter-to-iter
%               subspace angle history (Def: [])
%  FwRes      : Use Fast ADMM with Restart (Goldstein, et al. 2014)
%  FwResEta   : Fast ADMM with Restart parameter (Def: 0.999) 
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
%  by     Sejong Yoon (sjyoon@cs.rutgers.edu)
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
%      Alternating Direction Optimization Methods, SIAM J. Img. Sci. 7(3),
%      pp. 1588-1623, 2014.

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
defaultMaxIter = 1000;
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
addParameter(p,'ConvCounter', defaultConvCounter);
addParameter(p,'ConvCounterRMS', defaultConvCounterRMS);

parse(p,varargin{:});

% Initialize parameters
model_init  = p.Results.InitModel;
THRESH      = p.Results.Threshold;
THRESHA     = p.Results.ThresholdA;
THRESHR     = p.Results.ThresholdR;
iter_obj    = p.Results.ShowObjPer;
COUNTER_MAX = p.Results.MaxIter;
ZeroMean    = p.Results.ZeroMean;
ETA         = p.Results.Eta;
VaryETA     = p.Results.VaryETA;
VaryETA_mu  = p.Results.VaryETA_mu;
VaryETA_tau = p.Results.VaryETA_tau;
VaryETA_max = p.Results.VaryETA_max;
UseResidual = p.Results.UseResidual;
W_GT        = p.Results.W_GT;
FwRes       = p.Results.FwRes;
FwResEta    = p.Results.FwResEta;
converge_counter_max     = p.Results.ConvCounter;
converge_counter_rms_max = p.Results.ConvCounterRMS;

assert(ETA > 0, 'Learning rate (ETA) should be positive!');
assert(VaryETA_mu > 0, 'Please check 0 < VaryETA_mu < 1');
assert(VaryETA_mu < 1, 'Please check 0 < VaryETA_mu < 1');
assert(VaryETA_tau > 0, 'Please check 0 < VaryETA_tau');
assert(~(FwRes == true && VaryETA == true), ...
    'Please choose either VaryETA or FwRes, not both');

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
clearvars -global;

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

% dual variables
zbar_W = cell(J,1);
zbar_MU = cell(J,1);
zbar_PREC = cell(J,1);
o_zbar_W = cell(J,1);
o_zbar_MU = cell(J,1);
o_zbar_PREC = cell(J,1);

% residuals
pri_norms = cell(J,1);
dual_norms = cell(J,1);

% Learning rate
ETAhalf = ETA * 0.5;

% Initialize objective function - Lagrangian (we are minimizing this)
oldObjLR = realmax;
objArray = zeros(COUNTER_MAX, J+1); % last one is reserved for total
rArray = zeros(COUNTER_MAX, 1); 
sArray = zeros(COUNTER_MAX, 1); 
rtArray = zeros(COUNTER_MAX, 1); 
stArray = zeros(COUNTER_MAX, 1); 
ssaArray = zeros(COUNTER_MAX, 1);
rmsArray = zeros(COUNTER_MAX, 1);
alphaArray = zeros(COUNTER_MAX, J);
cArray = zeros(COUNTER_MAX, J);
Fi = zeros(J, 1);

%--------------------------------------------------------------------------
% Prepare performance measures
converged = 0;
counter = 1;
tic;

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
        Fi(i) = dppca_local( M, i, ETA, ZeroMean );
    end
    
    %----------------------------------------------------------------------
    % In each node: Update Lagrange multipliers
    for i = 1 : J
        Bi = Bj{i};
        for j = 1:length(Bi)
            LAMBDAi{i} = LAMBDAi{i} + ( ETAhalf * (Wi{i} - Wi{Bi(j)}) );
            GAMMAi{i} = GAMMAi{i} + ( ETAhalf * (MUi{i} - MUi{Bi(j)}) );
            BETAi{i} = BETAi{i} + ( ETAhalf * (PRECi{i} - PRECi{Bi(j)}) );
        end
    end
    
    %----------------------------------------------------------------------
    % Stopping criterion checkpoint

    % Compute objective
    objLR = 0;
    for i = 1 : J
        objLRi = Fi(i);
        Bi = Bj{i};
        for j = 1:length(Bi) 
            objLRi = objLRi ...
                + trace(LAMBDAi{i}' * (Wi{i} - Wi{Bi(j)})) ...
                + (GAMMAi{i}' * (MUi{i} - MUi{Bi(j)})) ...
                + (BETAi{i} * (PRECi{i} - PRECi{Bi(j)})) ...
                + ETAhalf * norm(Wi{i} - Wi{Bi(j)},'fro')^2 ...
                + ETAhalf * norm(MUi{i} - MUi{Bi(j)},'fro')^2 ...
                + ETAhalf * (PRECi{i} - PRECi{Bi(j)})^2;
        end
        objArray(counter, i) = objLRi;
        objLR = objLR + objLRi;
    end
    objArray(counter,J+1) = objLR;
    relErr = (objLR - oldObjLR) / abs(oldObjLR);
    oldObjLR = objLR;
    
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
        dual_tol_norm = sqrt(   norm(ETA * LAMBDAi{i},'fro')^2 ...
                              + norm(ETA * GAMMAi{i}, 'fro')^2 ...
                              + norm(ETA * BETAi{i},  'fro')^2     );
        dual_tol_i = THRESHA*sqrt(D*M + D + 1) + THRESHR*dual_tol_norm;
        dual_tol = dual_tol + dual_tol_i;

        % compute consensus dual residual
        dual_norm_i = sqrt(  norm(ETA * (zbar_W{i}    - o_zbar_W{i}),   'fro')^2 ...
                           + norm(ETA * (zbar_MU{i}   - o_zbar_MU{i}),  'fro')^2 ...
                           + norm(ETA * (zbar_PREC{i} - o_zbar_PREC{i}),'fro')^2    );
        dual_norm = dual_norm + dual_norm_i;
        
        % compute consensus primal residual
        pri_norm_i = sqrt(   norm(Wi{i}    - zbar_W{i},   'fro')^2 ...
                           + norm(MUi{i}   - zbar_MU{i},  'fro')^2 ...
                           + norm(PRECi{i} - zbar_PREC{i},'fro')^2     );
        pri_norm = pri_norm + pri_norm_i;

        % save
        pri_norms{i} = pri_norm_i;
        dual_norms{i} = dual_norm_i;
    end
    rArray(counter) = pri_norm;
    sArray(counter) = dual_norm;
    rtArray(counter) = pri_tol;
    stArray(counter) = dual_tol;
    
    % Adjust Eta (as He, et al. 2000)
    if VaryETA && counter <= VaryETA_max
        if pri_norm * VaryETA_mu >= dual_norm
            ETA = ETA * (1 + VaryETA_tau);
            ETAhalf = ETA * 0.5;
        elseif dual_norm * VaryETA_mu >= pri_norm
            ETA =  ETA / (1 + VaryETA_tau);
            ETAhalf = ETA * 0.5;
        end 
    end
    
    % Compute c_k
    if (FwRes)
        % compute c_k
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
                if iter_obj > 0 && mod(counter, iter_obj) == 0
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
                fprintf('  Node %d, c_{k} = %.2f, c_{k-1} = n/a, alpha_{k} = %.2f, alpha_{k-1} = n/a (INIT)\n', ...
                    i, cArray(counter,i), alphaArray(counter,i));
            end
        end
    end

    % Calculate subspace angle if requested
    if ~isempty(W_GT)
        ssaArray(counter) = calc_ppca_max_ssa_gt(Wi, W_GT, J);
    end    
    % Calculate reprojection error
    rmsArray(counter) = calc_ppca_rms(X, Wi, EZ, MUi);
    
    % Show progress if requested
    if iter_obj > 0 && mod(counter, iter_obj) == 0
        if UseResidual
            fprintf('Iter = %d:  p_norm = %f / p_tol = %f / d_norm = %f / d_tol = %f, RMS = %f, MSA = %.2e (J = %d, ETA = %f)\n', ...
                counter, pri_norm, pri_tol, dual_norm, dual_tol, ...
                calc_ppca_rms(Xj, Wi, EZ, MUi), ...
                calc_ppca_max_ssa(oWi, Wi), ...
                J, ETA);
        else
            fprintf('Iter = %d:  Cost = %f (rel %3.2f%%), RMS = %f, MSA = %.2e (J = %d, ETA = %f)\n', ...
                counter, objLR, relErr*100, ...
                calc_ppca_rms(Xj, Wi, EZ, MUi), ...
                calc_ppca_max_ssa(oWi, Wi), ...
                J, ETA);
        end
    end
    
    % Check whether it has converged
    if UseResidual
        if J == 1 && abs(relErr) < THRESH
            converged = 1;
            break;
        elseif J > 1 && pri_norm <= pri_tol && dual_norm <= dual_tol
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

    % Increase counter
    counter = counter + 1;
end

% Check convergence
if converged ~= 1
    counter = COUNTER_MAX;
    fprintf('Could not converge within %d iterations.\n', COUNTER_MAX);
end

% Compute performance measures
eTIME = toc;
eITER = counter;% Fill in remaining objective function value slots.
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

% Remove rotational ambiguity
%for idx = 1 : J
%     [~, ~, V] = svd(Wi(:,:,idx));
%     Wi(:,:,idx) = Wi(:,:,idx) * V;
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
%end

% Assign return values
W = cell(J,1);
MU = cell(J,1);
VAR = cell(J,1);
for i = 1 : J
    W{i} = Wi{i};
    MU{i} = MUi{i};
    VAR{i} = 1./PRECi{i};
end

% Create structure
model = structure( ...
    W, MU, VAR, ...
    EZ, EZZt, ...
    eITER, eTIME, objArray, rArray, sArray, rtArray, stArray, ...
    ssaArray, rmsArray, ...
    LAMBDAi, GAMMAi, BETAi);

% Clean up 
clearvars -global Xj;
clearvars -global Wi MUi PRECi oWi oMUi oPRECi;
clearvars -global EZ EZZt;
clearvars -global Bj MISSj;
clearvars -global LAMBDAi GAMMAi BETAi;

clearvars -except model;

end

function [F_new] = dppca_local( M, i, ETA, isMeanZero )
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

    W_new1 = sum( EZnZnt(:,:,NcI), 3 ) * oPRECi{i} + 2*ETA*cBj*eye(M);
    W_new2 = ( (Xi(d,NcI) - oMUi{i}(d)) * EZn(:,NcI)' ) * oPRECi{i};
    W_new3 = 2 * LAMBDAi{i}(d,:);
    W_new4 = zeros(1, M);
    for j = 1:cBj
        W_new4 = W_new4 + (oWi{i}(d,:) + oWi{Bi(j)}(d,:));
    end
    
    % Update
    if sum(sum(W_new1)) < eps
        W_new(d,:) = zeros(1,M);
    else
        W_new(d,:) = (W_new2 - W_new3 + ETA * W_new4) / W_new1;
    end
end

% Update MUi
MU_new = zeros(D,1);
if ~isMeanZero
    for d = 1 : D
        % Get non-missing point indexes
        NcI = (MISSi(d,:) == 0);
        
        MU_new1 = sum(NcI) * oPRECi{i} + 2*ETA*cBj;
        MU_new2 = oPRECi{i} * sum( Xi(d,NcI) - W_new(d,:) * EZn(:,NcI), 2 );
        MU_new3 = 2 * GAMMAi{i}(d);
        MU_new4 = 0;
        for j = 1 : cBj
            MU_new4 = MU_new4 + ( oMUi{i}(d) + oMUi{Bi(j)}(d) );
        end
        
        % Update
        MU_new(d) = (MU_new2 - MU_new3 + ETA * MU_new4) / MU_new1;
    end
end

% Update PRECi (by solve for VARi^(-1))
PREC_new1 = 2 * ETA * cBj;
PREC_new21 = 2 * BETAi{i};
PREC_new22 = 0;
for j = 1:cBj
    PREC_new22 = PREC_new22 + ETA * (oPRECi{i} + oPRECi{Bi(j)});
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

% Compute data log likelihood (we don't need to compute constant)
obj_val1 = 0;
obj_val2 = 0;
obj_val3 = 0;
obj_val4 = 0;
obj_val5 = 0;
for n = 1:Ni
    Id = (MISSi(:,n) == 0);
    Xc = Xi(Id,n);
    MUc = MU_new(Id);
    Wc = W_new(Id,:);        

    obj_val1 = obj_val1 + 0.5 * trace(EZnZnt(:,:,n));
    obj_val2 = obj_val2 + 0.5 * sum(Id) * log(2 * pi * PREC_new);
    obj_val3 = obj_val3 + (PREC_new/2) * norm(Xc - MUc, 'fro').^2;
    obj_val4 = obj_val4 + (PREC_new/2) * trace(EZnZnt(:,:,n) * (Wc' * Wc));
    obj_val5 = obj_val5 + PREC_new * EZn(:,n)' * Wc' * (Xc - MUc);
end
F_new = obj_val1 - obj_val2 + obj_val3 + obj_val4 - obj_val5;

% Update parameter values
Wi{i} = W_new;
MUi{i} = MU_new;
PRECi{i} = PREC_new;

EZ{i} = EZn;
EZZt{i} = EZnZnt;

end
