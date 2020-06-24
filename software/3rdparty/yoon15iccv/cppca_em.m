function model = cppca_em( X, M, varargin )
% CPPCA_EM   Probablistic Principal Component Analysis (PPCA) using EM
% 
% Description
%  Compute estimates of W and VAR given (D dimension)x(N observations) 
%  sample matrix X and convergence threshold using EM algorithm. Returns 
%  same output to cppca(X, M) but uses EM algorithm. We don't compute 
%  sample covariance matrix, eigenvalues and eigenvectors here. The code 
%  applied [3] to speed up convergence. NaN elements in X are considered 
%  as missing values.
%
% Input
%  X        : D x N observation matrix
%  M        : Scalar for latent space dimension
%
% Optional Parameters
%  InitModel  : PPCA model to set initial parameter (Def: random)
%  Threshold  : Scalar convergence criterion (Def: 1e-5)
%  ShowObjPer : If > 0, print out objective every specified iteration. 
%               If 0, nothing will be printed. (Def: 1)
%  MaxIter    : Maximum iterations (Def: 10000)
%  ZeroMean   : True if we enforce the mean to be zero. (Def: false)
%  W_GT       : Ground truth projection matrix if you want iter-to-iter
%               subspace angle history (Def: [])
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
%  on     2012.06.01 (last modified on 2015/02/13)
%
% References
%  [1] M.E. Tipping and C.M. Bishop, Probablistic principal component 
%      analysis, J. R. Statist. Soc. B (1999), 61 Part 3, pp. 611-622.
%  [2] Probablistic Modeling Toolkit 3, pmtk3.googlecode.com
%  [3] A. Ilin and T. Raiko, Practical Approaches to Principal Component
%      Analysis in the Presense of Missing VAlues, JMLR (2010), 1957-2000.

% Check required arguments
assert(nargin >= 2, 'Please specify at least X and M.');

% D x N = original dimension x samples
[D, N] = size(X);

% Set default convergence criterion
p = inputParser;
p.StructExpand = false;

W = orth(randn(D,M));
MU = zeros(D,1);
VAR = 1;
defaultMODEL = structure(W, MU, VAR);
defaultTHRESH = 1e-5;
defaultITER = 1;
defaultMaxIter = 10000;
defaultZeroMean = false;
defaultW_GT = [];

addParameter(p,'InitModel',defaultMODEL);
addParameter(p,'Threshold',defaultTHRESH,@isnumeric);
addParameter(p,'ShowObjPer',defaultITER,@isnumeric);
addParameter(p,'ZeroMean',defaultZeroMean);
addParameter(p,'MaxIter',defaultMaxIter);
addParameter(p,'W_GT',defaultW_GT);

parse(p,varargin{:});

% Initialize parameters
model_init  = p.Results.InitModel;
THRESH      = p.Results.Threshold;
iter_obj    = p.Results.ShowObjPer;
COUNTER_MAX = p.Results.MaxIter;
W_GT        = p.Results.W_GT;

% Get missing value indices
MissIDX = isnan(X);

%--------------------------------------------------------------------------
% Check validity of initilaization
if (isfield(model_init, 'W') && iscell(model_init.W)) || ...
    (isfield(model_init, 'MU') && iscell(model_init.MU)) || ...
    (isfield(model_init, 'VAR') && iscell(model_init.VAR))
    error('Invalid initialization: please specify centralized model');
end

% Initialize latent variables (for loop implementation)
EZn = randn(M, N);
EZnZnt = repmat(eye(M, M), [1, 1, N]);

% Initialize parameters; Note that we don't need to initialize MU
W = model_init.W;
MU = model_init.MU; 
VAR = model_init.VAR;

% Initialize data log likelihood as 0 and prepare the constant
oldLL = -realmax;
objArray = zeros(COUNTER_MAX,1);
ssaArray = zeros(COUNTER_MAX,1);

% Prepare performance measures
converged = 0;
counter = 1;
tic;

% Main loop
while counter <= COUNTER_MAX
    %----------------------------------------------------------------------
    % E step
    %----------------------------------------------------------------------

    for n = 1:N    
        % Get indices of available features
        Id = (MissIDX(:,n) == 0);
    
        % Compute Minv = (W'W + VAR*I) first
        Wc = W(Id,:);
        Minv = Wc'*Wc + VAR*eye(M);
    
        % E[z_n] = Minv * W' * (x_n - MU)
        % Currently M x N
        EZn(:,n) = Minv \ Wc' * (X(Id,n) - MU(Id));

        % E[z_n z_n'] = VAR * Minv + E[z_n]E[z_n]'
        % Currently M x M
        EZnZnt(:,:,n) = VAR * inv(Minv) + EZn(:,n) * EZn(:,n)';
    end
    
    %----------------------------------------------------------------------
    % M step
    %----------------------------------------------------------------------
    
    % W_new (Eq. 12.56, PRML (Bishop, 2006) pp.578)
    W_new = zeros(D,M);
    for d = 1:D
        % Get indices of samples with current feature
        In = (MissIDX(d,:) == 0);
        if sum(In) == 0, continue; end
        
        W_new1 = (X(d,In) - MU(d)) * EZn(:,In)';
        W_new2 = sum(EZnZnt(:,:,In), 3);
        W_new(d,:) = W_new1 /  W_new2;
    end
    
    % MU_new (NOTE: this is only being assigned here to make the recursion 
    % identical to that of the distributed PPCA)
    MU_new = zeros(D,1);
    if ~p.Results.ZeroMean
        for d = 1:D
            % Get indices of samples with current feature
            In = (MissIDX(d,:) == 0);
            if sum(In) == 0, continue; end

            MU_new(d) = mean( X(d,In) - W_new(d,:)*EZn(:,In), 2 );
        end    
    end

    % VAR_new (Eq. 12.57, PRML (Bishop, 2006) pp.578)
    VAR_new1 = 0;
    VAR_new2 = 0;
    VAR_new3 = 0;
    VAR_new4 = 0;
    for n = 1:N
        % Get indices of available features
        Id = (MissIDX(:,n) == 0);
        VAR_new1 = VAR_new1 + sum(Id);
        
        VAR_new2 = VAR_new2 + norm((X(Id,n) - MU_new(Id)), 2)^2;
        VAR_new3 = VAR_new3 + 2*(EZn(:,n)' * W_new(Id,:)' * (X(Id,n) - MU_new(Id)));
        VAR_new4 = VAR_new4 + trace(EZnZnt(:,:,n) * (W_new(Id,:)' * W_new(Id,:)));
    end
    VAR_new1 = 1 / VAR_new1;
    VAR_new = VAR_new1 * (VAR_new2 - VAR_new3 + VAR_new4);

    % Correct rotational ambiguity (Sec. 4.1, 4.2, Ilin and Raiko, 2010)
    % To check Eq. 34 hold, try:
    %   mean(model.EZ,2)
    % To check Eq. 35 hold, try:
    %   Sig = zeros(2,2); 
    %   for n = 1:N, 
    %     Sig = Sig + (model.EZ(:,n)*model.EZ(:,n)' + model.EZZt(:,:,n)); 
    %   end; 
    %   Sig = Sig ./ N; Sig
    mEZn = mean(EZn,2);
    MU_new = MU_new + W_new * mEZn;
    EZn = bsxfun(@minus, EZn, mEZn);
    
    vEZnZnt = EZn * EZn';
    for n = 1:N
        vEZnZnt = vEZnZnt + EZnZnt(:,:,n);
    end
    vEZnZnt = vEZnZnt ./ N;
    [Uz, Dz] = eig(vEZnZnt);
    Dz = sqrt(Dz);
    vWWt = ((W_new * Uz * Dz)' * (W_new * Uz * Dz)) ./ D;
    [Vw, Dw] = eig(vWWt);
    [~, I] = sort(-diag(Dw));
    Vw = Vw(:,I);

    W_new = W_new * Uz * Dz * Vw;
    R = Vw' * diag(1 ./ diag(Dz)) * Uz';
    EZn = R * EZn;
    for n = 1:N
        EZnZnt(:,:,n) = R * EZnZnt(:,:,n) * R';
    end
    
    %----------------------------------------------------------------------
    % Check convergence
    %----------------------------------------------------------------------
    
    % Compute data log likelihood (we don't need to compute constant)
    obj_val1 = 0;
    obj_val2 = 0;
    obj_val3 = 0;
    obj_val4 = 0;
    obj_val5 = 0;
    for n = 1:N
        Id = (MissIDX(:,n) == 0);
        Xc = X(Id,n);
        MUc = MU_new(Id);
        Wc = W_new(Id,:);        
        
        obj_val1 = obj_val1 + 0.5 * sum(Id) * log(2 * pi * VAR_new);
        obj_val2 = obj_val2 + 0.5 * trace(EZnZnt(:,:,n));
        obj_val3 = obj_val3 + (1/(2*VAR_new)) * norm(Xc - MUc, 2).^2;
        obj_val4 = obj_val4 + (1/VAR_new) * EZn(:,n)' * Wc' * (Xc - MUc);
        obj_val5 = obj_val5 + (1/(2*VAR_new)) * trace(EZnZnt(:,:,n) * (Wc' * Wc));
    end
    LL = -(obj_val1 + obj_val2 + obj_val3 - obj_val4 + obj_val5);
    
    objArray(counter,1) = LL;
    relErr = (LL - oldLL)/abs(oldLL);
    oldLL = LL;
    
    % Calculate subspace angle if requested
    if ~isempty(W_GT)
        ssaArray(counter) = calc_ppca_max_ssa_gt(W_new, W_GT);
    end
    
    % Show progress if requested
    if iter_obj > 0 && (mod(counter, iter_obj) == 0)
        fprintf('Iter %d: LL = %f (rel %3.2f%%), RMS = %f, SA = %.2e\n', ...
            counter, ...
            LL, relErr*100, ...
            calc_ppca_rms(X, W_new, EZn, MU_new), ...
            calc_ppca_max_ssa(W, W_new)); 
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

% Truncate to save memory
ssaArray = ssaArray(1:counter);

% Remove rotational ambiguity (we already did this within the main loop)
% It seems like this is not enough; we need to enforce these conditions: 
%   Z ~ N(0,I) and W' * W = diag(s) where s_k >= s_l if k < l.
% [~, ~, V] = svd(W);
% W = W * V;
% EZ = V'*EZ;
% for idn=1:N
%     EZZt(:,:,idn) = V'*EZZt(:,:,idn)*V;
% end

% Remove rotational ambiguity for SfM
% [W, EZ, Qinv] = remove_rot_amb(W, EZ);
% for idn = 1 : N
%     EZZt(:,:,idn) = Qinv * EZZt(:,:,idn) * Qinv';
% end

% % Compute optimization formula (why we need log probability)
% S = bsxfun(@minus, X, MU);
% S = S * S' / N;
% C = W * W' + VAR * eye(D);
% F_new = -(N/2) * (D*log(2*pi) + log(det(C)) + trace(C\S));

% Create model structure
model = structure(W, MU, VAR, EZ, EZZt, eITER, eTIME, objArray, ssaArray);

end
