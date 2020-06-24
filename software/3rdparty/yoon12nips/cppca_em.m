function model = cppca_em( X, M, varargin )
% CPPCA_EM   Probablistic Principal Component Analysis (PPCA) using EM
% 
% Description
%  Compute estimates of W and VAR given (D dimension)x(N observations) 
%  sample matrix X and convergence threshold using EM algorithm. Returns 
%  same output to cppca(X, M) but uses EM algorithm thus W*W' will be the
%  equal not W. We don't compute sample covariance matrix, eigenvalues and 
%  eigenvectors here. We used loops deliverately for clearer presentation 
%  by sacrificing speed. 
%
% Input
%  X          : D x N observation matrix
%  M          : Scalar for latent space dimension
%  [Optional Parameters]
%  InitModel  : PPCA model to set initial parameter (Def: random)
%  Threshold  : Scalar convergence criterion (Def: 1e-6)
%  ShowObjPer : If > 0, print out objective every specified iteration. 
%               If 0, nothing will be printed. (Def: 10)
%  ZeroMean   : True if we enforce the mean to be zero. (Def: false)
%  MissIDX    : D x N binary matrix indicating a feature is missing, i.e.
%               true means the feature is missing. (Def: all zeros)
%
% Output
%  model = structure(W, MU, VAR, EZ, EZZt, eITER, eTIME, objArray)
%  W        : D x M projection matrix
%  MU       : D x 1 vector sample means
%  VAR      : Scalar estimated variance
%  EZ       : M x N matrix mean of N latent vectors
%  EZZt     : M x M x N cube for covariance of N latent vectors
%  eITER    : Iterations took
%  eTIME    : Elapsed time
%  objArray : Objective function value change over iterations
%
% Implemented
%  by     Sejong Yoon (sjyoon@cs.rutgers.edu)
%  on     2012.06.01 (last modified on 2014/06/05)
%
% References
%  [1] M.E. Tipping and C.M. Bishop, Probablistic principal component 
%      analysis, J. R. Statist. Soc. B (1999), 61 Part 3, pp. 611-622.
%  [2] Probablistic Modeling Toolkit 3, pmtk3.googlecode.com

COUNTER_MAX = 10000;

% Check required arguments
assert(nargin >= 2, 'Please specify at least X and M.');

% D x N = original dimension x samples
[D, N] = size(X);

% Set default convergence criterion
p = inputParser;
p.StructExpand = false;

W = rand(D,M);
MU = rand(D,1);
VAR = rand(1);
defaultMODEL = structure(W, MU, VAR);
defaultTHRESH = 1e-6;
defaultITER = 10;
defaultZeroMean = false;
defaultMissIDX = zeros(D,N);

addParamValue(p,'InitModel',defaultMODEL);
addParamValue(p,'Threshold',defaultTHRESH,@isnumeric);
addParamValue(p,'ShowObjPer',defaultITER,@isnumeric);
addParamValue(p,'ZeroMean',defaultZeroMean);
addParamValue(p,'MissIDX',defaultMissIDX);

parse(p,varargin{:});

% Initialize parameters
% Note that we don't need to initialize MU
W = p.Results.InitModel.W;
MU = p.Results.InitModel.MU; 
VAR = p.Results.InitModel.VAR;

THRESH = p.Results.Threshold;
iter_obj = p.Results.ShowObjPer;
MissIDX = p.Results.MissIDX;

% Initialize latent variables (for loop implementation)
EZn = zeros(M,N);
EZnZnt = zeros(M,M,N);

% Initialize data log likelihood as 0 and prepare the constant
oldLL = -realmax;
objArray = zeros(COUNTER_MAX,1);

% Prepare performance measures
converged = 0;
counter = 1;
tic;

% Main loop
while counter <= COUNTER_MAX
    %----------------------------------------------------------------------
    % E step
    %----------------------------------------------------------------------

    for idn = 1:N    
        % Get indices of available features
        DcI = (MissIDX(:,idn) == 0);
    
        % Compute Minv = (W'W + VAR*I) first
        Wc = W(DcI,:);
        Minv = Wc'*Wc + VAR*eye(M);
    
        % E[z_n] = Minv * W' * (x_n - MU)
        % Currently M x N
        EZn(:,idn) = Minv \ Wc' * (X(DcI,idn) - MU(DcI));

        % E[z_n z_n'] = VAR * Minv + E[z_n]E[z_n]'
        % Currently M x M
        EZnZnt(:,:,idn) = VAR * inv(Minv) + EZn(:,idn) * EZn(:,idn)';
    end
    
    %----------------------------------------------------------------------
    % M step
    %----------------------------------------------------------------------
    
    % W_new (Eq. 12.56, PRML (Bishop, 2006) pp.578)
    W_new = zeros(D,M);
    for idd = 1:D
        % Get indices of samples with current feature
        NcI = (MissIDX(idd,:) == 0);
        if sum(NcI) == 0, continue; end
        
        W_new1 = (X(idd,NcI) - MU(idd)) * EZn(:,NcI)';
        W_new2 = sum(EZnZnt(:,:,NcI), 3);
        W_new(idd,:) = W_new1 /  W_new2;
    end
    
    % MU_new  (NOTE: this is only being assigned here to make the recursion 
    % identical to that of the distributed PCCA)
    MU_new = zeros(D,1);
    if ~p.Results.ZeroMean
        for idd = 1:D
            % Get indices of samples with current feature
            NcI = (MissIDX(idd,:) == 0);
            if sum(NcI) == 0, continue; end

            MU_new(idd) = mean( X(idd,NcI) - W_new(idd,:)*EZn(:,NcI), 2 );
        end    
    end

    % VAR_new (Eq. 12.57, PRML (Bishop, 2006) pp.578)
    VAR_new1 = 0;
    VAR_new2 = 0;
    VAR_new3 = 0;
    VAR_new4 = 0;
    for idn = 1:N
        % Get indices of available features
        DcI = (MissIDX(:,idn) == 0);
        VAR_new1 = VAR_new1 + sum(DcI);
        
        VAR_new2 = VAR_new2 + norm((X(DcI,idn) - MU_new(DcI)), 2)^2;
        VAR_new3 = VAR_new3 + 2*(EZn(:,idn)' * W_new(DcI,:)' * (X(DcI,idn) - MU_new(DcI)));
        VAR_new4 = VAR_new4 + trace(EZnZnt(:,:,idn) * (W_new(DcI,:)' * W_new(DcI,:)));
    end
    VAR_new1 = 1 / VAR_new1;
    VAR_new = VAR_new1 * (VAR_new2 - VAR_new3 + VAR_new4);
    
    %----------------------------------------------------------------------
    % Check convergence
    %----------------------------------------------------------------------
    
    % Compute data log likelihood (we don't need to compute constant)
    obj_val1 = 0;
    obj_val2 = 0;
    obj_val3 = 0;
    obj_val4 = 0;
    obj_val5 = 0;
    for idn = 1:N
        DcI = (MissIDX(:,idn) == 0);
        Xc = X(DcI,idn);
        MUc = MU_new(DcI);
        Wc = W_new(DcI,:);        
        
        obj_val1 = obj_val1 + 0.5 * sum(DcI) * log(2 * pi * VAR_new);
        obj_val2 = obj_val2 + 0.5 * trace(EZnZnt(:,:,idn));
        obj_val3 = obj_val3 + (1/(2*VAR_new)) * norm(Xc - MUc, 2).^2;
        obj_val4 = obj_val4 + (1/VAR_new) * EZn(:,idn)' * Wc' * (Xc - MUc);
        obj_val5 = obj_val5 + (1/(2*VAR_new)) * trace(EZnZnt(:,:,idn) * (Wc' * Wc));
    end
    LL = -(obj_val1 + obj_val2 + obj_val3 - obj_val4 + obj_val5);
    
    objArray(counter,1) = LL;
    relErr = (LL - oldLL)/abs(oldLL);
    oldLL = LL;
    
    % Show progress if requested
    if nargin >= 5 && iter_obj > 0 && (mod(counter, iter_obj) == 0)
        fprintf('Iter %d:  LL = %f (rel %3.2f%%)\n', counter, LL, relErr*100); 
    end

    % Update parameters and latent statistics with new values
    W = W_new;
    VAR = VAR_new;
    MU = MU_new;
    EZ = EZn;
    EZZt = EZnZnt;
    
    % Check convergence
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
eITER = counter;
eTIME = toc;

% Fill in remaining objective function value slots. Convert objective to
% minimization function so as to compare with distributed version
if counter < COUNTER_MAX
    objArray(counter+1:COUNTER_MAX,1) = ...
        repmat(objArray(counter,1), [COUNTER_MAX - counter, 1]);
end
objArray = objArray .* -1;

% Remove rotational ambiguity
[~, ~, V] = svd(W);
W = W * V;
EZ = V'*EZ;
for idn=1:N
    EZZt(:,:,idn) = V'*EZZt(:,:,idn)*V;
end
% Remove rotational ambiguity for SfM
% [W, EZ, Qinv] = remove_rot_amb(W, EZ);
% for idn = 1 : N
%     EZZt(:,:,idn) = Qinv * EZZt(:,:,idn) * Qinv';
% end

% Create model structure
model = structure(W, MU, VAR, EZ, EZZt, eITER, eTIME, objArray);

end
