function model = dppca( X, M, V, E, varargin )
% DPPCA      Distributed Probablistic PCA (D-PPCA)
% 
% Description
%  Solve probabilistic PCA problem in a distributed way. The network has 
%  max(V) nodes. We assume the network is connected. This function only 
%  simulates parameter broadcasts. Local computation is done by dppca_local 
%  function.
%
% Input
%  X     : D x N matrix for full data from all nodes (N=N1+N2+...+NJ)
%  M     : Scalar of projection dimension
%  V     : N x 1 vector for each observation's source (node affiliation)
%  E     : J x J adjacency matrix where J = max(V)
%  [Optional Parameters]
%  InitModel  : D-PPCA model to set initial parameter (Def: random)
%  Threshold  : Scalar convergence criterion (Def: 10^-5)
%  ShowObjPer : If > 0, print out objective every specified iteration.
%               If 0, nothing will be printed. (Def: 10)
%  ZeroMean   : True if we enforce the mean to be zero. (Def: false)
%  MissIDX    : D x N binary matrix indicating a feature is missing, i.e.
%               true means the feature is missing. (Def: all zeros)
%  Eta        : Scalar of learning rate (Def: 10)
%
% Output
%  model = structure(W, MU, VAR, ...);
%  W        : D x M x J projection matrices for J nodes
%  MU       : D x J vector sample means for J nodes
%  VAR      : J x 1 scalar estimated variances for J nodes
%  EZ       : J cells; M x N matrix, mean of N latent vectors
%  EZZt     : J cells; M x M x N cube, covariance of N latent vectors
%  eITER    : Iterations took
%  eTIME    : Elapsed time
%  objArray : Objective function value change over iterations
%  LAMBDAi, GAMMAi, BETAi : Lagrange multipliers for debugging purpose
%
% Implemented
%  by     Sejong Yoon (sjyoon@cs.rutgers.edu)
%  on     2011.10.07 (last modified on 2014/08/19)
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

COUNTER_MAX = 10000;

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

W = rand(D, M, J);
MU = rand(D, 1);
VAR = rand(1)./10;
defaultMODEL = structure(W, MU, VAR);
defaultTHRESH = 1e-5;
defaultITER = 10;
defaultZeroMean = false;
defaultMissIDX = zeros(D,N);

defaultETA = 10;

addParamValue(p,'InitModel',defaultMODEL);
addParamValue(p,'Threshold',defaultTHRESH,@isnumeric);
addParamValue(p,'ShowObjPer',defaultITER,@isnumeric);
addParamValue(p,'ZeroMean',defaultZeroMean);
addParamValue(p,'MissIDX',defaultMissIDX);

addParamValue(p,'Eta',defaultETA,@isnumeric);

parse(p,varargin{:});

% Initialize parameters
model_init = p.Results.InitModel;
THRESH     = p.Results.Threshold;
iter_obj   = p.Results.ShowObjPer;
ZeroMean   = p.Results.ZeroMean;
MissIDX    = p.Results.MissIDX;

ETA        = p.Results.Eta;

assert(ETA > 0, 'Learning rate (ETA) should be positive!');

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
for idx = 1 : J
    Xj{idx} = X(:, V == idx);
end

% Find i-th node's neighbor set Bi in advance to speed up
Bj = cell(J,1);
for idx = 1 : J
    if r == 0
        Bj{idx} = [];
    else
        Bj{idx} = find(E(idx,:) > 0);
    end
end

% Local parameters and auxiliary variables defined here for simplicity.
% In the real environment, these variables reside in each local sensor.
% Initialize parameters
Wi = zeros(D, M, J);
for idx = 1 : J
    Wi(:,:,idx) = model_init.W(:,:,idx);
end
MUi = repmat(model_init.MU, [1, J]);
PRECi = repmat(1./model_init.VAR, [1, J]);

% Initialize Lagrange multipliers. Each edge of each node has a multiplier.
LAMBDAi = zeros(D, M, J);
GAMMAi = zeros(D, J);
BETAi = zeros(J, 1);

% Initialize latent variables
EZ = cell(J, 1);
EZZt = cell(J, 1);

% Build MISSi for speed up
MISSj = cell(J,1);
for idx = 1 : J
    MISSj{idx} = MissIDX(:, V == idx);
end

% Learning rate
ETAhalf = ETA * 0.5;

% Initialize objective function - Lagrangian (we are minimizing this)
oldObjLR = realmax;
objArray = zeros(COUNTER_MAX, J+1); % last one is reserved for total
Fi = zeros(J, 1);

%--------------------------------------------------------------------------
% Prepare performance measures
converged = 0;
counter = 1;
tic;

% Main loop
while counter <= COUNTER_MAX
    %----------------------------------------------------------------------
    % Temporarily store parameters to simulate broadcasting and
    % synchronization. All nodes should update before sending their values.
    oWi = Wi;
    oMUi = MUi;
    oPRECi = PRECi;
    
    %----------------------------------------------------------------------
    % In each node: Update parameters locally
    for idx = 1 : J
        Fi(idx) = dppca_local( M, idx, ETA, ZeroMean );
    end
    
    %----------------------------------------------------------------------
    % In each node: Update Lagrange multipliers
    for idx = 1 : J
        Bi = Bj{idx};
        for jn = 1:length(Bi)
            LAMBDAi(:,:,idx) = LAMBDAi(:,:,idx) + ...
                ( ETAhalf * (Wi(:,:,idx) - Wi(:,:,Bi(jn))) );
            GAMMAi(:,idx) = GAMMAi(:,idx) + ...
                ( ETAhalf * (MUi(:,idx) - MUi(:,Bi(jn))) );
            BETAi(idx) = BETAi(idx) + ...
                ( ETAhalf * (PRECi(idx) - PRECi(Bi(jn))) );
        end
    end
    
    %----------------------------------------------------------------------
    % Stopping criterion checkpoint
        
    % Compute objective
    objLR = 0;
    for idx = 1 : J
        objLRi = Fi(idx);
        Bi = Bj{idx};
        for jn = 1:length(Bi) 
            objLRi = objLRi ...
                + trace(LAMBDAi(:,:,idx)' * ...
                    (Wi(:,:,idx) - Wi(:,:,Bi(jn)))) ...
                + (GAMMAi(:,idx)' * ...
                    (MUi(:,idx) - MUi(:,Bi(jn)))) ...
                + (BETAi(idx) * ...
                    (PRECi(idx) - PRECi(Bi(jn)))) ...
                + ETAhalf * ...
                    norm(Wi(:,:,idx) - Wi(:,:,Bi(jn)),2)^2 ...
                + ETAhalf * ...
                    norm(MUi(:,idx) - MUi(:,Bi(jn)),2)^2 ...
                + ETAhalf * ...
                    (PRECi(idx)-PRECi(Bi(jn)))^2;
        end
        objArray(counter, idx) = objLRi;
        objLR = objLR + objLRi;
    end
    objArray(counter,J+1) = objLR;
    relErr = (objLR - oldObjLR) / abs(oldObjLR);
    oldObjLR = objLR;
    
    % Show progress if requested
    if iter_obj > 0 && mod(counter, iter_obj) == 0
        fprintf('Iter = %d:  Cost = %f (rel %3.2f%%) (J = %d, ETA = %f)\n', ...
            counter, objLR, relErr*100, J, ETA);
    end
    
    % Check whether it has converged
    if abs(relErr) < THRESH
        converged = 1;
        break;
    end

    % Increase counter
    counter = counter + 1;
end

% Check convergence
if converged ~= 1
    fprintf('Could not converge within %d iterations.\n', COUNTER_MAX);
end

% Compute performance measures
eTIME = toc;
eITER = counter;

% Fill in remaining objective function value slots.
if counter < COUNTER_MAX
    % fill in the remaining objArray
    for idx = 1:J+1
        objArray(counter+1:COUNTER_MAX,idx) = ...
            repmat(objArray(counter,idx), [COUNTER_MAX - counter, 1]);
    end
end

% Remove rotational ambiguity
for idx = 1 : J
    [~, ~, V] = svd(Wi(:,:,idx));
    Wi(:,:,idx) = Wi(:,:,idx) * V;
    EZ{idx} = V'*EZ{idx};
    for idy = 1 : size(Xj{idx},2)
        EZZt{idx}(:,:,idy) = V' * EZZt{idx}(:,:,idy) * V;
    end
%     % OR...
%     [t1, t2] = remove_rot_amb(EZ{idx}', Wi(:,:,idx)');
%     EZ{idx} = t1';
%     Wi(:,:,idx) = t2';
end


% Assign return values
W = Wi;
MU = MUi;
VAR = 1./PRECi;

% Create structure
model = structure( ...
    W, MU, VAR, ...
    EZ, EZZt, ...
    eITER, eTIME, objArray, ...
    LAMBDAi, GAMMAi, BETAi);

clearvars -except model;

end

function [F_new] = dppca_local( M, idx, ETA, isMeanZero )
% DPPCA_LOCAL  D-PPCA Local Update
% 
% Input
%  M     : Projected dimension
%  idx   : Current node index
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
Xi = Xj{idx};
Bi = Bj{idx};
MISSi = MISSj{idx};

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
    Wc = Wi(DcI,:,idx);
    Mi = Wc' * Wc + 1/PRECi(idx) * eye(M);

    % E[Zn] = Mi^(-1) * Wi' * (Xin - MUi)
    % Currently M x N
    EZn(:,n) = Mi \ Wc' * (Xi(DcI,n) - MUi(DcI,idx));

    % E[z_n z_n'] = VAR * Mi^(-1) + E[z_n]E[z_n]'
    % Currently M x M
    EZnZnt(:,:,n) = inv(Mi) * (1/PRECi(idx)) + EZn(:,n) * EZn(:,n)';
end

%--------------------------------------------------------------------------
% M-step

% Update Wi
W_new = zeros(D, M);
% One dimension at a time
for d = 1 : D
    % Get non-missing point indexes
    NcI = (MISSi(d,:) == 0);

    W_new1 = sum( EZnZnt(:,:,NcI), 3 ) * oPRECi(idx) + 2*ETA*cBj*eye(M);
    W_new2 = ( (Xi(d,NcI) - oMUi(d,idx)) * EZn(:,NcI)' ) * oPRECi(idx);
    W_new3 = 2 * LAMBDAi(d,:,idx);
    W_new4 = zeros(1, M);
    for jn = 1:cBj
        W_new4 = W_new4 + (oWi(d,:,idx) + oWi(d,:,Bi(jn)));
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
        
        MU_new1 = sum(NcI) * oPRECi(idx) + 2*ETA*cBj;
        MU_new2 = oPRECi(idx) * sum( Xi(d,NcI) - W_new(d,:) * EZn(:,NcI), 2 );
        MU_new3 = 2 * GAMMAi(d,idx);
        MU_new4 = 0;
        for jn = 1 : cBj
            MU_new4 = MU_new4 + ( oMUi(d,idx) + oMUi(d,Bi(jn)) );
        end
        
        % Update
        MU_new(d) = (MU_new2 - MU_new3 + ETA * MU_new4) / MU_new1;
    end
end

% Update PRECi (by solve for VARi^(-1))
PREC_new1 = 2 * ETA * cBj;
PREC_new21 = 2 * BETAi(idx);
PREC_new22 = 0;
for jn = 1:cBj
    PREC_new22 = PREC_new22 + ETA * (oPRECi(idx) + oPRECi(Bi(jn)));
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
        + trace( EZnZnt(:,:,n) * Wc' * Wc ) );
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

% Compute optimization formula
S = bsxfun(@minus, Xi, MU_new);
S = S * S' / Ni;
C = W_new * W_new' + (1/PREC_new) * eye(D);
if D > M
    logDetC = (D-M)*log(1/PREC_new)+log(det(W_new'*W_new+(1/PREC_new)*eye(M)));
else
    logDetC = log(det(C));
end
F_new = (Ni/2)*(logDetC + trace(C\S));

% Update parameter values
Wi(:,:,idx) = W_new;
MUi(:,idx) = MU_new;
PRECi(idx) = PREC_new;

EZ{idx} = EZn;
EZZt{idx} = EZnZnt;

end
