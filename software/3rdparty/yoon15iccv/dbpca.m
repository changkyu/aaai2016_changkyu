function model = dbpca( X, M, V, E, varargin )
% DBPCA      Distributed Variational Bayesian PCA (D-BPCA)
% 
% Description
%  Solve Variational Bayes PCA problem in a distributed way. The network 
%  has max(V) nodes. We assume the network is connected. This function only 
%  simulates parameter broadcasts. Local computation is done by dbpca_local 
%  function. NaN elements in X are considered as missing values.
%
% Input
%  X     : D x N matrix for full data from all nodes (N=N1+N2+...+NJ)
%  M     : Scalar of projection dimension
%  V     : N x 1 vector for each observation's source (node affiliation)
%  E     : J x J adjacency matrix where J = max(V)
%  [Optional Parameters]
%  InitModel  : D-BPCA model to set initial parameter (Def: random)
%  Threshold  : Scalar convergence criterion (Def: 1e-5)
%  ShowObjPer : If > 0, print out objective every specified iteration.
%               If 0, nothing will be printed. (Def: 1)
%  MaxIter    : Maximum iterations (Def: 10000)
%  ZeroMean   : True if we enforce the mean to be zero. (Def: false)
%  Eta        : Scalar of learning rate (Def: 10)
%
% Output
%  model = structure(...);
%  mW  : J cells; D x M matrix, mean     of W
%  vW  : J cells; D x M matrix, variance of W
%  mMU : J cells; D x 1 vector, mean     of MU
%  vMU : J cells; D x 1 matrix, variance of MU
%  PX  : J cells; Noise precision parameters (data)   / scalar
%  PW  : J cells; Noise precision parameters (weight) / D x M matrix
%  PMU : J cells; Noise precision parameters (mean)   / scalar
%  mZ  : J cells; M x N matrix, mean of N latent vectors
%  vZ  : J cells; M x M x N cube, variance of N latent vectors
%  eITER    : Iterations took
%  eTIME    : Elapsed time
%  objArray : Objective function value change over iterations
%  gW  : J cells; D x M matrix, Lagrange multiplier for mW
%  bW  : J cells; D x M matrix, Lagrange multiplier for vW
%  gMU : J cells; D x 1 vector, Lagrange multiplier for mMU
%  bMU : J cells; D x 1 vector, Lagrange multiplier for vMU
%  gPX : J cells; Lagrange multiplier for PX  / scalar
%  gPW : J cells; Lagrange multiplier for PW  / D x M matrix
%  gPMU: J cells; Lagrange multiplier for PMU / scalar
%
% Implemented
%  by     Sejong Yoon (sjyoon@cs.rutgers.edu)
%  on     2014.12.01 (last modified on 2015/03/24)
%
% References
%  [1] TBA

% Check required arguments
assert(nargin >= 2, 'Please specify at least X and M.');

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
% We need to broadcast parameters & Lagrange multipliers to neighbors. 
% Here, we use global variables. In real settings, sensors should transmit 
% them over network.

% Data and Model parameters
global Xj;
global mZ vZ mW vW mMU vMU;
global PX PW PMU;
global omW ovW omMU ovMU;
global oPX oPW oPMU;
global barW barMU;
global Bj MISSj;
global gMU bMU gW bW gPX gPW gPMU;

%--------------------------------------------------------------------------
% Parse optional parameters
p = inputParser;
p.StructExpand = false;

mW = cell(J, 1);
mMU = cell(J, 1);
PX = cell(J, 1);
PW = cell(J, 1);
PMU = cell(J, 1);
barW = cell(J, 1);
barMU = cell(J, 1);
for i = 1 : J
    mW{i} = orth(randn(D, M));
    mMU{i} = zeros(D, 1);
    PX{i} = 1;
    PW{i} = ones(D,M);
    PMU{i} = 1;
    barW{i} = zeros(D, M);
    barMU{i} = zeros(D, 1);
end
defaultMODEL = structure(mW, mMU, PX, PW, PMU, barW, barMU);
defaultTHRESH = 1e-5;
defaultITER = 1;
defaultMaxIter = 10000;
defaultZeroMean = false;
defaultETA = 10;

addParameter(p,'InitModel',defaultMODEL);
addParameter(p,'Threshold',defaultTHRESH,@isnumeric);
addParameter(p,'ShowObjPer',defaultITER,@isnumeric);
addParameter(p,'MaxIter',defaultMaxIter);
addParameter(p,'ZeroMean',defaultZeroMean);
addParameter(p,'Eta',defaultETA,@isnumeric);

parse(p,varargin{:});

% Initialize parameters
model_init = p.Results.InitModel;
THRESH     = p.Results.Threshold;
iter_obj   = p.Results.ShowObjPer;
COUNTER_MAX = p.Results.MaxIter;
ZeroMean   = p.Results.ZeroMean;
ETA        = p.Results.Eta;

assert(ETA > 0, 'Learning rate (ETA) should be positive!');

%--------------------------------------------------------------------------
% Check validity of initilaization
if (isfield(model_init, 'W') && ~iscell(model_init.W)) || ...
    (isfield(model_init, 'MU') && ~iscell(model_init.MU)) || ...
    (isfield(model_init, 'VAR') && ~iscell(model_init.VAR)) || ...        
    (isfield(model_init, 'mW') && ~iscell(model_init.mW)) || ...
    (isfield(model_init, 'mMU') && ~iscell(model_init.mMU)) || ...
    (isfield(model_init, 'PX') && ~iscell(model_init.PX)) || ...
    (isfield(model_init, 'PW') && ~iscell(model_init.PW)) || ...
    (isfield(model_init, 'PMU') && ~iscell(model_init.PMU)) || ...
    (isfield(model_init, 'barW') && ~iscell(model_init.barW)) || ...
    (isfield(model_init, 'barMU') && ~iscell(model_init.barMU)),
    error('Invalid initialization: please specify distributed model');
end

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

% Initialize latent variables and model parameters
mZ = cell(J, 1);
vZ = cell(J, 1);
for i = 1 : J
    mZ{i} = randn(M, N);
    vZ{i} = repmat(eye(M, M), [1, 1, N]);
end

if isfield(model_init, 'mW')
    mW = model_init.mW;
elseif isfield(model_init, 'W')
    mW = model_init.W;
else
    error('Invalid initialization: need either W or mW');
end
for i = 1 : J
    vW{i} = ones(D, M);
end

if isfield(model_init, 'mMU')
    mMU = model_init.mMU;
elseif isfield(model_init, 'MU')
    mMU = model_init.MU;
else
    error('Invalid initialization: need either MU or mMU');
end
for i = 1 : J
    vMU{i} = ones(D, 1);
end

if isfield(model_init, 'PX') && isfield(model_init, 'PW') && isfield(model_init, 'PMU')
    PX = model_init.PX;
    PW = repmat(model_init.PW, [D, M]);
    PMU = model_init.PMU;
elseif isfield(model_init, 'VAR')
    PX = model_init.VAR;
    PW = cell(J, 1);
    PMU = cell(J, 1);
    for i = 1 : J
        PW{i} = ones(D, M);
        PMU{i} = 1;
    end
else
    error('Invalid initialization: need either VAR or PX');
end

if isfield(model_init, 'barW')
    barW = model_init.barW;
end
if isfield(model_init, 'barMU')
    barMU = model_init.barMU;
end

% Build MISSi for speed up
MISSj = cell(J,1);
for i = 1 : J
    MISSj{i} = isnan(Xj{i});
end

% Initialize Lagrange multipliers. Each edge of each node has a multiplier.
gW = cell(J, 1);
bW = cell(J, 1);
gMU = cell(J, 1);
bMU = cell(J, 1);
gPX = cell(J, 1);
gPW = cell(J, 1);
gPMU = cell(J, 1);
for i = 1 : J
    gW{i} = zeros(D, M);
    bW{i} = zeros(D, M);
    gMU{i} = zeros(D, 1);
    bMU{i} = zeros(D, 1);
    gPX{i} = 0;
    gPW{i} = zeros(D, M);
    gPMU{i} = 0;
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
    omW = mW;
    ovW = vW;
    omMU = mMU;
    ovMU = vMU;
    oPX = PX;
    oPW = PW;
    oPMU = PMU;
    
    %----------------------------------------------------------------------
    % In each node: Update parameters locally
    for i = 1 : J
        Fi(i) = dbpca_local( M, i, ETA, ZeroMean );
    end
        
    %----------------------------------------------------------------------
    % In each node: Update Lagrange multipliers
    for i = 1 : J
        Bi = Bj{i};
        for j = 1 : length(Bi)
            % Eq. 86
            gW{i} = gW{i} + ETAhalf * ( mW{i} - mW{Bi(j)} );
            % Eq. 61
            gMU{i} = gMU{i} + ETAhalf * ( mMU{i} - mMU{Bi(j)} );
            % Eq. 64
            bW{i} = bW{i} + ETAhalf * ( vW{i} - vW{Bi(j)} );
            % Eq. 62
            bMU{i} = bMU{i} + ETAhalf * ( vMU{i} - vMU{Bi(j)} );
            % Eq. 65
            gPX{i} = gPX{i} + ETAhalf * ( PX{i} - PX{Bi(j)} );
            % Eq. 66
            gPW{i} = gPW{i} + ETAhalf * ( PW{i} - PW{Bi(j)} );
            % Eq. 67
            gPMU{i} = gPMU{i} + ETAhalf * ( PMU{i} - PMU{Bi(j)} );
        end
    end
    
    %----------------------------------------------------------------------
    % Stopping criterion checkpoint
    
    % Compute objective (Eq. 10)
    objLR = 0;
    for i = 1 : J
        objLRi = Fi(i);
        Bi = Bj{i};
        for j = 1 : length(Bi)
            % W-related
            for d = 1 : D
                for m = 1 : M
                    objLRi = objLRi ...
                        + gW{i}(d,m) * (mW{i}(d,m) - mW{Bi(j)}(d,m)) ...
                        + bW{i}(d,m) * (vW{i}(d,m) - vW{Bi(j)}(d,m)) ...
                        + gPW{i}(d,m) * (PW{i}(d,m) - PW{Bi(j)}(d,m)) ...
                        + ETAhalf * (mW{i}(d,m) - mW{Bi(j)}(d,m)).^2 ...
                        + ETAhalf * (vW{i}(d,m) - vW{Bi(j)}(d,m)).^2 ...
                        + ETAhalf * (PW{i}(d,m) - PW{Bi(j)}(d,m)).^2;
                end
            end
            % Others
            objLRi = objLRi ...
                + gMU{i}' * (mMU{i} - mMU{Bi(j)}) ...
                + bMU{i}' * (vMU{i} - vMU{Bi(j)}) ...
                + gPX{i} * (PX{i} - PX{Bi(j)}) ...
                + gPMU{i} * (PMU{i} - PMU{Bi(j)}) ...
                + ETAhalf * norm(mMU{i} - mMU{Bi(j)},2).^2 ...
                + ETAhalf * norm(vMU{i} - vMU{Bi(j)},2).^2 ...
                + ETAhalf * (PX{i} - PX{Bi(j)}).^2 ...
                + ETAhalf * (PMU{i} - PMU{Bi(j)}).^2;
        end
        objArray(counter, i) = objLRi;
        objLR = objLR + objLRi;
    end
    objArray(counter, i+1) = objLR;
    relErr = (objLR - oldObjLR) / abs(oldObjLR);
    oldObjLR = objLR;
    
    % Show progress if requested
    if iter_obj > 0 && mod(counter, iter_obj) == 0
        fprintf('Iter %d: Cost = %f (rel %3.2f%%), RMS = %f, MSA = %.2e (J = %d, ETA = %f)\n', ...
            counter, ...
            objLR, relErr*100, ...
            calc_ppca_rms(Xj, mW, mZ, mMU), ...
            calc_ppca_max_ssa(omW, mW), ...
            J, ETA);
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
    for i = 1:J+1
        objArray(counter+1:COUNTER_MAX,i) = ...
            repmat(objArray(counter,i), [COUNTER_MAX - counter, 1]);
    end
end

% Convert precision to variance
for i = 1 : J
    PX{i} = 1./PX{i};
    PW{i} = 1./PW{i};
    PMU{i} = 1./PMU{i};

    % Invert vZ from precision to variance
    Ni = size(vZ{i},3);
    for n = 1 : Ni
        vZ{i}(:,:,n) = inv(vZ{i}(:,:,n));
    end
end

% Create structure
model = structure( ...
    PX, PW, PMU, ...
    mZ, vZ, mW, vW, mMU, vMU, ...
    barW, barMU, ...
    eITER, eTIME, objArray, ...
    gMU, bMU, gW, bW, gPX, gPW, gPMU);

% Clean up 
clearvars -global Xj;
clearvars -global mZ vZ mW vW mMU vMU;
clearvars -global PX PW PMU;
clearvars -global omW ovW omMU ovMU;
clearvars -global oPX oPW oPMU;
clearvars -global barW barMU;
clearvars -global Bj MISSj;
clearvars -global gMU bMU gW bW gPX gPW gPMU;

clearvars -except model;

end

%%
function [F_new] = dbpca_local( M, i, ETA, ZeroMean )
% DBPCA_LOCAL  BPCA Local Update
% 
% Input
%  M        : Projected dimension
%  i        : Current node index
%  ETA      : Learning rate
%  ZeroMean : True if we fix mean
%
% Output
%  F_new    : 1 x 1 scalar computed optimization forumla (first term only)

% Parameters and latent variables. Other variables are defined as global 
% just for simple and easy-to-understand implementation.
global Xj;
global mZ vZ mW vW mMU vMU;
global PX PW PMU;
global omW ovW omMU ovMU;
global oPX oPW oPMU;
global barW barMU;
global Bj MISSj;
global gMU bMU gW bW gPX gPW gPMU;

% Take i-th node
Xi = Xj{i};
Bi = Bj{i};
MISSi = MISSj{i};

% Get size of this samples and ball of this node
[D, Ni] = size(Xi);
cBj = length(Bi);

% Initialize latent variables (for loop implementation)
mZ{i} = zeros(M, Ni);
vZ{i} = zeros(M, M, Ni);

%--------------------------------------------------------------------------
% E-step (eq. 13, 14)
%--------------------------------------------------------------------------
for n = 1 : Ni
    % Get indices of available features for this sample
    Id = (MISSi(:,n) == 0);
    
    % eq. 13
    vZ{i}(:,:,n) = eye(M, M) + oPX{i} * (omW{i}' * omW{i} + diag(sum(ovW{i},1)));

    % eq. 14
    mZ{i}(:,n) = oPX{i} * (vZ{i}(:,:,n) \ (omW{i}(Id,:)' * (Xi(Id,n) - omMU{i}(Id))));
end

%--------------------------------------------------------------------------
% M-step
%--------------------------------------------------------------------------
% Update vMU (eq. 76)
MUcoeff0 = -0.5;
for d = 1 : D
    MUcoeff1 = 0.5 * (Ni * oPX{i} + oPMU{i}) + 2 * bMU{i}(d);
    for j = 1 : cBj
        MUcoeff1 = MUcoeff1 - ETA * (ovMU{i}(d) + ovMU{Bi(j)}(d));
    end
    MUcoeff2 = 2 * ETA * cBj;
    
    MUnew = roots([MUcoeff2, MUcoeff1, MUcoeff0]);
    vMU{i}(d) = get_best_from_solve(MUnew);
end

% Update mMU (eq. 81)
mMU{i} = zeros(D, 1);
if ~ZeroMean
    MUdenom = (Ni * oPX{i} + oPMU{i} + 2 * ETA * cBj) * eye(D);

    MUnum = zeros(D, 1);
    for d = 1 : D
        In = (MISSi(d,:) == 0);

        MUnum(d) = MUnum(d) + ...
            sum( oPX{i} * ( Xi(d,In) - omW{i}(d,:) * mZ{i}(:,In) ) );
    end
    MUnum = MUnum - 2 * gMU{i};
    for j = 1 : cBj
        MUnum = MUnum + ETA * (omMU{i} + omMU{Bi(j)});
    end
    MUnum = MUnum + oPMU{i} * barMU{i};

    mMU{i} = MUdenom \ MUnum;
end

% Update vW (eq.79)
Pi = zeros(M, M);
for n = 1 : Ni
    Pi = Pi + ( mZ{i}(:,n) * mZ{i}(:,n)' + inv(vZ{i}(:,:,n)) );
end

Wcoeff0 = -0.5;
for d = 1 : D
    for m = 1 : M
        Wcoeff1 = (oPX{i} * Pi(m,m) + oPW{i}(d,m)) / 2 + 2 * bW{i}(d,m);
        for j = 1 : cBj
            Wcoeff1 = Wcoeff1 - ETA * (ovW{i}(d,m) + ovW{Bi(j)}(d,m));
        end
        Wcoeff2 = 2 * ETA * cBj;

        Wnew = roots([Wcoeff2, Wcoeff1, Wcoeff0]);
        vW{i}(d,m) = get_best_from_solve(Wnew);
    end
end

% Update mW (eq.85)
Wnum = zeros(D, M);
for n = 1 : Ni
    Id = (MISSi(:,n) == 0);
    
    Wnum(Id,:) = Wnum(Id,:) + (Xi(Id,n) - mMU{i}(Id)) * mZ{i}(:,n)';
end
Wnum = oPX{i} * Wnum;

for d = 1 : D
    Wnum(d,:) = Wnum(d,:) - 2 * gW{i}(d,:) + oPW{i}(d,:) .* barW{i}(d,:);
    for j = 1 : cBj
        Wnum(d,:) = Wnum(d,:) +  ETA * (omW{i}(d,:) + omW{Bi(j)}(d,:));
    end

    Wdenom = zeros(M, M);
    for n = 1 : Ni
        Wdenom = Wdenom + (mZ{i}(:,n) * mZ{i}(:,n)' + inv(vZ{i}(:,:,n)));
    end
    Wdenom = oPX{i} * Wdenom;
    Wdenom = Wdenom + diag(oPW{i}(d,:)) + (2 * ETA * cBj) * eye(M);

    mW{i}(d,:) = Wnum(d,:) / Wdenom;
end

% Compute A(t+1) (eq.71)
A = 0;
for n = 1 : Ni
    % Get indices of available features for this sample
    Id = (MISSi(:,n) == 0);
    
    A = A + (Xi(Id,n)' * Xi(Id,n));
    A = A - 2 * (Xi(Id,n)' * mW{i}(Id,:) * mZ{i}(:,n));
    A = A - 2 * (Xi(Id,n)' * mMU{i}(Id));
    A = A + 2 * (mMU{i}' * mW{i} * mZ{i}(:,n));
    A = A + (mMU{i}' * mMU{i} + sum(vMU{i},1));
    A = A + trace( ( mZ{i}(:,n) * mZ{i}(:,n)' + inv(vZ{i}(:,:,n)) ) ...
        * ( mW{i}' * mW{i} + diag(sum(vW{i},1)) ) );
end

% Update PX (eq.72)
PXcoeff0 = -Ni * D / 2;
PXcoeff1 = A / 2 + 2 * gPX{i};
for j = 1 : cBj
    PXcoeff1 = PXcoeff1 - ETA * (oPX{i} + oPX{Bi(j)});
end
PXcoeff2 = 2 * ETA * cBj;

Pnew = roots([PXcoeff2, PXcoeff1, PXcoeff0]);
PX{i} = get_best_from_solve( Pnew );

% Update PW (eq.73)
for d = 1 : D
    for m = 1 : M
        PWcoeff0 = -1/2;
        PWcoeff1 = 0.5 * ( mW{i}(d,m).^2 + vW{i}(d,m) + barW{i}(d,m).^2 ...
            - 2 * mW{i}(d,m) * barW{i}(d,m) ) + 2 * gPW{i}(d,m);
        for j = 1 : cBj
            PWcoeff1 = PWcoeff1 - ETA * (oPW{i}(d,m) + oPW{Bi(j)}(d,m));
        end
        PWcoeff2 = 2 * ETA * cBj;

        Pnew = roots([PWcoeff2, PWcoeff1, PWcoeff0]);
        PW{i}(d,m) = get_best_from_solve( Pnew );
    end
end

% Update PMU (eq.74)
PMUcoeff0 = -D / 2;
PMUcoeff1 = 0.5 * ( (mMU{i}' * mMU{i}) + sum(vMU{i},1) ...
    + (barMU{i}' * barMU{i}) - 2 * (mMU{i}' * barMU{i}) ) + 2 * gPMU{i};
for j = 1 : cBj
    PMUcoeff1 = PMUcoeff1 - ETA * (oPMU{i} + oPMU{Bi(j)});
end
PMUcoeff2 = 2 * ETA * cBj;

Pnew = roots([PMUcoeff2, PMUcoeff1, PMUcoeff0]);
PMU{i} = get_best_from_solve(Pnew);

%--------------------------------------------------------------------------
% Correct rotational ambiguity (Sec. 4.1, 4.2, Ilin and Raiko, 2010)
% mEZ = mean(mZ{i},2);
% mMU{i} = mMU{i} + mW{i} * mEZ;
% mZ{i} = bsxfun(@minus, mZ{i}, mEZ);
% 
% vEZnZnt = mZ{i} * mZ{i}';
% for idn=1:Ni
%     vEZnZnt = vEZnZnt + inv(vZ{i}(:,:,idn)); % vZ is precision here
% end
% vEZnZnt = vEZnZnt ./ Ni;
% [Uz, Dz] = eig(vEZnZnt);
% Dz = sqrt(Dz);
% vWWt = ((mW{i} * Uz * Dz)' * (mW{i} * Uz * Dz));
% for idd=1:D
%     vW{i}(idd,:) = diag((Uz * Dz)' * diag(vW{i}(idd,:)) * (Uz * Dz));
%     vWWt = vWWt + diag(vW{i}(idd,:));
% end
% vWWt = vWWt ./ D;
% [Vw, Dw] = eig(vWWt);
% [~, I] = sort(-diag(Dw));
% Vw = Vw(:,I);
% 
% mW{i} = mW{i} * Uz * Dz * Vw;
% for idd=1:D
%     vW{i}(idd,:) = diag(Vw' * diag(vW{i}(idd,:)) * Vw);
% end
% R = Vw' * diag(1 ./ diag(Dz)) * Uz';
% mZ{i} = R * mZ{i};
% for idn=1:Ni
%     vZ{i}(:,:,idn) = inv((R / vZ{i}(:,:,idn)) * R'); % vZ is precision here
% end

%--------------------------------------------------------------------------
% Compute optimization formula (eq. 10)
t01 = A; % actually this is the same as A we computed above
% t01 = 0;
% for n = 1 : Ni
%     t01 = t01 + Xi(:,n)' * Xi(:,n);
%     t01 = t01 - 2 * Xi(:,n)' * mW{idx} * mZ{idx}(:,n);
%     t01 = t01 - 2 * Xi(:,n)' * mMU{idx};
%     t01 = t01 + 2 * mMU{idx}' * mW{idx} * mZ{idx}(:,n);
%     t01 = t01 + mMU{idx}' * mMU{idx} + sum(vMU{idx},1);
%     t01 = t01 + trace( ( mZ{idx}(:,n) * mZ{idx}(:,n)' + inv(vZ{idx}(:,:,n)) ) ...
%         * ( mW{idx}' * mW{idx} + diag(sum(vW{idx},1)) ) );
% end
t01 = t01 * PX{i} / 2;

t02 = -(Ni * D / 2) * log(PX{i});

t03 = 0;
for n = 1 : Ni
    t03 = t03 + ( mZ{i}(:,n)' * mZ{i}(:,n) + sum(diag(inv(vZ{i}(:,:,n)))) );
end
t03 = t03 / 2;

t04 = 0;
for d = 1 : D
    for m = 1 : M
        t04 = t04 + PW{i}(d,m) ...
            * ( (mW{i}(d,m) - barW{i}(d,m)) * (mW{i}(d,m) - barW{i}(d,m))' + vW{i}(d,m) );
    end
end
t04 = t04 / 2;

t05 = -(1 / 2) * sum(sum(log(PW{i})));
t06 = (PMU{i} / 2) * ( (mMU{i} - barMU{i})' * (mMU{i} - barMU{i}) + sum(vMU{i},1) );
t07 = -(D / 2) * log(PMU{i});
t08 = 0;
for n = 1 : Ni
    t08 = t08 + log( det( vZ{i}(:,:,n) ) );
end
t08 = t08 / 2;
t09 = -0.5 * sum(sum(log(vW{i})));
t10 = -0.5 * sum(log(vMU{i}));

F_new = t01 + t02 + t03 + t04 + t05 + t06 + t07 + t08 + t09 + t10;

end

%% Helper function to choose (hopefully) the best solution (positive reals)
function sol = get_best_from_solve(sols)

if length(sols) == 1
    sol = sols(1);
elseif length(sols) == 2
    if isreal(sols(1))
        if isreal(sols(2))
            sol = max(sols);
        else
            sol = sols(1);
        end
    elseif isreal(sols(2))
        sol = sols(2);
    else
        error('No real solution found!');
    end
else
    sol = -Inf;
    for r = 1 : length(sols)
        if isreal(sols(r)) && (sol < sols(r))
            sol = sols(r);
        end
    end
end

end
