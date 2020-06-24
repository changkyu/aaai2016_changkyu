function model = cbpca( X, M, varargin )
% CBPCA      Centralized Variational Bayesian PCA (C-BPCA)
% 
% Description
%  Solve Bayesian PCA problem in using VB method. NaN elements in X are 
%  considered as missing values.
%
% Input
%  X     : D x N matrix for full data
%  M     : Scalar of projection dimension
%  [Optional Parameters]
%  InitModel  : C-BPCA model to set initial parameter (Def: random)
%  Threshold  : Scalar convergence criterion (Def: 1e-5)
%  ShowObjPer : If > 0, print out objective every specified iteration.
%               If 0, nothing will be printed. (Def: 1)
%  MaxIter    : Maximum iterations (Def: 10000)
%  ZeroMean   : True if we enforce the mean to be zero. (Def: false)
%
% Output
%  model = structure(...);
%  PX  : Noise precision parameters (data)   / scalar
%  PW  : Noise precision parameters (weight) / D x M matrix
%  PMU : Noise precision parameters (mean)   / scalar
%  mZ  : M x N matrix, mean of N latent vectors
%  vZ  : M x M x N cube, variance of N latent vectors
%  mW  : D x M matrix, mean of W
%  vW  : D x M matrix, variance of W
%  mMU : D x 1 vector, mean of Mu
%  vMU : D x 1 matrix, variance of Mu
%  eITER    : Iterations took
%  eTIME    : Elapsed time
%  objArray : Objective function value change over iterations
%
% Implemented
%  by     Sejong Yoon (sjyoon@cs.rutgers.edu)
%  on     2014.12.01 (last modified on 2014/12/04)
%
% References
%  [1] TBA

% Check required arguments
assert(nargin >= 2, 'Please specify at least X and M.');

% D dimensions x N samples
[D, N] = size(X);

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
global MissIDX;

%--------------------------------------------------------------------------
% Parse optional parameters
p = inputParser;
p.StructExpand = false;

mW = orth(rand(D, M));
mMU = zeros(D, 1);
PX = 1;
PW = ones(D, M);
PMU = 1;
barW = zeros(D, M);
barMU = zeros(D, 1);
defaultMODEL = structure(mW, mMU, PX, PW, PMU, barW, barMU);
defaultTHRESH = 1e-5;
defaultITER = 1;
defaultMaxIter = 10000;
defaultZeroMean = false;

addParameter(p,'InitModel',defaultMODEL);
addParameter(p,'Threshold',defaultTHRESH,@isnumeric);
addParameter(p,'ShowObjPer',defaultITER,@isnumeric);
addParameter(p,'MaxIter',defaultMaxIter);
addParameter(p,'ZeroMean',defaultZeroMean);

parse(p,varargin{:});

% Initialize parameters
model_init  = p.Results.InitModel;
THRESH      = p.Results.Threshold;
iter_obj    = p.Results.ShowObjPer;
COUNTER_MAX = p.Results.MaxIter;
ZeroMean    = p.Results.ZeroMean;

% Get missing value indices
MissIDX = isnan(X);

%--------------------------------------------------------------------------
% Check validity of initilaization
if (isfield(model_init, 'W') && iscell(model_init.W)) || ...
    (isfield(model_init, 'MU') && iscell(model_init.MU)) || ...
    (isfield(model_init, 'VAR') && iscell(model_init.VAR)) || ...        
    (isfield(model_init, 'mW') && iscell(model_init.mW)) || ...
    (isfield(model_init, 'mMU') && iscell(model_init.mMU)) || ...
    (isfield(model_init, 'PX') && iscell(model_init.PX)) || ...
    (isfield(model_init, 'PW') && iscell(model_init.PW)) || ...
    (isfield(model_init, 'PMU') && iscell(model_init.PMU))
    error('Invalid initialization: please specify centralized model');
end

% Build Xi for speed up
Xj = X;

% Initialize latent variables and model parameters
mZ = randn(M, N);
vZ = repmat(eye(M, M), [1, 1, N]);

if isfield(model_init, 'mW')
    mW = model_init.mW;
elseif isfield(model_init, 'W')
    mW = model_init.W;
else
    error('Invalid initialization: need either W or mW');
end
vW = ones(D, M);

if isfield(model_init, 'mMU')
    mMU = model_init.mMU;
elseif isfield(model_init, 'MU')
    mMU = model_init.MU;
else
    error('Invalid initialization: need either MU or mMU');
end
vMU = ones(D, 1);

if isfield(model_init, 'PX') && isfield(model_init, 'PW') && isfield(model_init, 'PMU')
    PX = model_init.PX;
    PW = model_init.PW;
    PMU = model_init.PMU;
elseif isfield(model_init, 'VAR')
    PX = model_init.VAR;
    PW = ones(D, M);
    PMU = 1;
else
    error('Invalid initialization: need either VAR or PX');
end

if isfield(model_init, 'barW')
    barW = model_init.barW;
else
    barW = zeros(size(mW));
end
if isfield(model_init, 'barMU')
    barMU = model_init.barMU;
else
    barMU = zeros(size(mMU));
end

% Initialize objective function - Lagrangian (we are minimizing this)
oldObjLR = realmax;
objArray = zeros(COUNTER_MAX, 1); % last one is reserved for total

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
    Fi = cbpca_local( M, ZeroMean );
        
    % Compute objective
    objLR = Fi;
    objArray(counter) = objLR;
    relErr = (objLR - oldObjLR) / abs(oldObjLR);
    oldObjLR = objLR;
    
    % Show progress if requested
    if iter_obj > 0 && mod(counter, iter_obj) == 0
        fprintf('Iter %d: Cost = %f (rel %3.2f%%), RMS = %f, SA = %.2e\n', ...
            counter, ...
            objLR, relErr*100, ...
            calc_ppca_rms(Xj, mW, mZ, mMU), ...
            calc_ppca_max_ssa(omW, mW));
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

% Correct variance
%[vW, vMU, vZ] = mfvb_correct(Xj, mW, mZ, mMU, vW, vZ, vMU, PX, PW, PMU);

% Compute performance measures
eTIME = toc;
eITER = counter;

% Fill in remaining objective function value slots.
if counter < COUNTER_MAX
    % fill in the remaining objArray
    objArray(counter+1:COUNTER_MAX) = ...
        repmat(objArray(counter), [COUNTER_MAX - counter, 1]);
end

% Convert precision to variance
PX = 1./PX;
PW = 1./PW;
PMU = 1./PMU;

% Invert vZ from precision to variance
for n = 1 : N
    vZ(:,:,n) = inv(vZ(:,:,n));
end

% Create structure
model = structure( ...
    PX, PW, PMU, ...
    mZ, vZ, mW, vW, mMU, vMU, ...
    barW, barMU, ...
    eITER, eTIME, objArray);

clearvars -except model;

% Clean up 
clearvars -global Xj;
clearvars -global mZ vZ mW vW mMU vMU;
clearvars -global PX PW PMU;
clearvars -global omW ovW omMU ovMU;
clearvars -global oPX oPW oPMU;
clearvars -global barW barMU;
clearvars -global MissIDX;

end

%%
function [F_new] = cbpca_local( M, ZeroMean )
% CBPCA_LOCAL  BPCA Local Update
% 
% Input
%  M        : Projected dimension
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
global MissIDX;

% Take i-th node
Xi = Xj;

% Get size of this samples and ball of this node
[D, N] = size(Xi);

% Initialize latent variables (for loop implementation)
mZ = zeros(M, N);
vZ = zeros(M, M, N);

%--------------------------------------------------------------------------
% E-step (eq. 13, 14)
%--------------------------------------------------------------------------
for n = 1 : N
    % Get indices of available features for this sample
    Id = (MissIDX(:,n) == 0);
    
    % eq. 13
    vZ(:,:,n) = eye(M, M) + oPX * (omW' * omW + diag(sum(ovW,1)));

    % eq. 14
    mZ(:,n) = oPX * (vZ(:,:,n) \ (omW(Id,:)' * (Xi(Id,n) - omMU(Id))));
end

%--------------------------------------------------------------------------
% M-step
%--------------------------------------------------------------------------
% Update vMU (eq. 76)
MUcoeff0 = -0.5;
for d = 1 : D
    MUcoeff1 = 0.5 * (N * oPX + oPMU);
    MUnew = roots([MUcoeff1, MUcoeff0]);
    vMU(d) = MUnew;
end

% Update mMU (eq. 81)
mMU = zeros(D, 1);
if ~ZeroMean
    MUdenom = (N * oPX + oPMU);
    MUdenom = MUdenom * eye(D);

    MUnum = zeros(D, 1);
    for d = 1 : D
        In = (MissIDX(d,:) == 0);

        MUnum(d) = MUnum(d) + sum( oPX * ( Xi(d,In) - mW(d,:) * mZ(:,In) ) );
    end
    MUnum = MUnum + oPMU * barMU;

    mMU = MUdenom \ MUnum;
end

% Update vW (eq.79)
Wcoeff0 = -0.5;

Pi = zeros(M, M);
for n = 1 : N
    Pi = Pi + ( mZ(:,n) * mZ(:,n)' + inv(vZ(:,:,n)) );
end

for d = 1 : D
    for m = 1 : M
        Wcoeff1 = (oPX * Pi(m,m) + oPW(d,m)) / 2;
        Wnew = roots([Wcoeff1, Wcoeff0]);
        vW(d,m) = Wnew;
    end
end

% Update mW (eq.84)
term1 = zeros(D, M);
for n = 1 : N
    Id = (MissIDX(:,n) == 0);
    
    term1(Id,:) = term1(Id,:) + (Xi(Id,n) - mMU(Id)) * mZ(:,n)';
end

for d = 1 : D
    term1t = oPX * term1(d,:) + oPW(d,:) .* barW(d,:);

    term2 = zeros(M, M);
    for n = 1 : N
        term2 = term2 + (mZ(:,n) * mZ(:,n)' + inv(vZ(:,:,n)));
    end
    term2 = oPX * term2;
    term2 = term2 + diag(oPW(d,:));

    mW(d,:) = term1t / term2;
end

% Compute A(t+1) (eq.71)
A = 0;
for n = 1 : N
    % Get indices of available features for this sample
    Id = (MissIDX(:,n) == 0);
    
    A = A + (Xi(Id,n)' * Xi(Id,n));
    A = A - 2 * (Xi(Id,n)' * mW(Id,:) * mZ(:,n));
    A = A - 2 * (Xi(Id,n)' * mMU(Id));
    A = A + 2 * (mMU' * mW * mZ(:,n));
    A = A + (mMU' * mMU + sum(vMU,1));
    A = A + trace( ( mZ(:,n) * mZ(:,n)' + inv(vZ(:,:,n)) ) ...
        * ( mW' * mW + diag(sum(vW,1)) ) );
end

% Update PX (eq.72)
PXcoeff0 = -N * D / 2;
PXcoeff1 = A / 2;

Pnew = roots([PXcoeff1, PXcoeff0]);
PX = Pnew;

% Update PW (eq.73)
for d = 1 : D
    for m = 1 : M
        PWcoeff0 = -1/2;
        PWcoeff1 = 0.5 * ( mW(d,m).^2 + vW(d,m) + barW(d,m).^2 - 2 * mW(d,m) * barW(d,m) );

        Pnew = roots([PWcoeff1, PWcoeff0]);
        PW(d,m) = max(Pnew);
    end
end

% Update PMU (eq.74)
PMUcoeff0 = -D / 2;
PMUcoeff1 = 0.5 * ( (mMU' * mMU) + sum(vMU,1) + (barMU' * barMU) - 2 * (mMU' * barMU) );

Pnew = roots([PMUcoeff1, PMUcoeff0]);
PMU = max(Pnew);

%--------------------------------------------------------------------------
% Correct rotational ambiguity (Sec. 4.1, 4.2, Ilin and Raiko, JMLR 2010)
mEZ = mean(mZ,2);
mMU = mMU + mW * mEZ;
mZ = bsxfun(@minus, mZ, mEZ);

vEZnZnt = mZ * mZ';
for idn=1:N
    vEZnZnt = vEZnZnt + inv(vZ(:,:,idn)); % vZ is precision here
end
vEZnZnt = vEZnZnt ./ N;
[Uz, Dz] = eig(vEZnZnt);
Dz = sqrt(Dz);
vWWt = ((mW * Uz * Dz)' * (mW * Uz * Dz));
for idd=1:D
    vW(idd,:) = diag((Uz * Dz)' * diag(vW(idd,:)) * (Uz * Dz));
    vWWt = vWWt + diag(vW(idd,:));
end
vWWt = vWWt ./ D;
[Vw, Dw] = eig(vWWt);
[~, I] = sort(-diag(Dw));
Vw = Vw(:,I);

mW = mW * Uz * Dz * Vw;
for idd=1:D
    vW(idd,:) = diag(Vw' * diag(vW(idd,:)) * Vw);
end
R = Vw' * diag(1 ./ diag(Dz)) * Uz';
mZ = R * mZ;
for idn=1:N
    vZ(:,:,idn) = inv((R / vZ(:,:,idn)) * R'); % vZ is precision here
end

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
t01 = t01 * PX / 2;

t02 = -(N * D / 2) * log(PX);

t03 = 0;
for n = 1 : N
    t03 = t03 + ( mZ(:,n)' * mZ(:,n) + sum(diag(inv(vZ(:,:,n)))) );
end
t03 = t03 / 2;

t04 = 0;
for d = 1 : D
    for m = 1 : M
        t04 = t04 + PW(d,m) * ( (mW(d,m) - barW(d,m)) * (mW(d,m) - barW(d,m))' + vW(d,m) );
    end
end
t04 = t04 / 2;

t05 = -(1 / 2) * sum(sum(log(PW)));
t06 = (PMU / 2) * ( (mMU - barMU)' * (mMU - barMU) + sum(vMU,1) );
t07 = -(D / 2) * log(PMU);
t08 = 0;
for n = 1 : N
    t08 = t08 + log( det( vZ(:,:,n) ) );
end
t08 = t08 / 2;
t09 = -0.5 * sum(sum(log(vW)));
t10 = -0.5 * sum(log(vMU));

F_new = t01 + t02 + t03 + t04 + t05 + t06 + t07 + t08 + t09 + t10;

end
