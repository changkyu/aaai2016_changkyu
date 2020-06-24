%  PCA_FULL - PCA with unrestricted Gaussian posterior
%
%  PCA with unrestricted Gaussian pdfs (full covariance matrices) for
%  approximating the posterior distributions of A(i,:) and S(:,j) in
%  the model X(:,j) = Mu + A*S(:,j) + Noise. The noise is isotropic
%  and Gaussian with the variance V.
%
%  [ A, S, Mu, V, CV, HP, LC ] = PCA_FULL( X, N ) identifies the model
%  for the given data matrix X and number of principal components N.
%  Matrix X can be either sparse with only observed values included
%  (note that observed zeros should be replaced with eps) or a full
%  matrix with missing values replaced by NaNs.
%
%  The function returns the mean values of the model parameters in A,
%  S, Mu, and V. The posterior variances are stored in CV such that
%  CV.A{j} is the posterior covariance matrix for A(i,:) and
%  CV.S{CV.Isv(j)} (or CV.S{j} if CV.Isv is empty) is the same for
%  S(:,j). HP contains point estimates of the hyperparameters: HP.Va
%  is the prior variance for A(:,k) and HP.Vmu is the same for Mu.
%  
%  LC is a structure with learning curves (rms training and probing
%  errors, cost function values, time).
%
%  PCA_FULL( X, N, 'algorithm', name ) specifies the algorithm:
%   'ppca': Probabilistic PCA (no prior and point estimates for A, Mu)
%   'map':  Prior on A(:,k) and Mu, point estimates for A(i,:) and Mu
%   'vb':   Variational Bayesian PCA (prior on A and Mu, Gaussian
%           posterior approximation for A(i,:) and Mu) {default}
%  
%  Other optional parameter/value pairs with {default values}:
%   init       - Initialization type ({'random'} - random,
%                filename: load from file, structure: from given data)
%   maxiters   - Maximum number of iterations {1000}
%   rotate2pca - Whether to perform rotation of A and S to speed-up
%                convergence {1}
%   uniquesv   - Whether to compute only unique covariance matrices of
%                S(:,j) {1}
%   minangle   - Termination by minimum angle between subspaces
%                defined by A {1e-8}
%   rmsstop    - Termination by rms training error ([] - no rms stop)
%                {[ 100 1e-4 1e-3 ]}: 100 iterations, 1e-4 absolute
%                tolerance, 1e-3 relative tolerance,
%   cfstop     - Termination by cost function value {[]} (similarly to
%                rmsstop). The cost function is computed only if this
%                option is nonempty
%   xprobe     - Validation data set (of the same dimensions as X)
%   earlystop  - Whether to use early stopping based on probing error
%   verbose    - Progress information display level (0,{1},2)
%   display    - Plot progress {0}
%   autosave   - Auto-save after each {600} seconds
%   filename   - Name of the file for auto-save {'pca_f_autosave'}
%
%  OUT = PCA_FULL( X, N ) returns all the outputs in a single
%  structure variable OUT. Learning can be continued as follows:
%    out = pca_full( X, N );
%    [ A, S, Mu, V, CV, HP, LC ] = pca_full( X, N, 'init', out );

%  TODO: More efficient implementation of V

%  This software is provided "as is", without warranty of any kind.
%  Alexander Ilin, Tapani Raiko

function [ A, S, Mu, V, cv, hp, lc ] = pca_full( X, ncomp, varargin )

opts = struct( ...
    'init',          'random',...
    'maxiters',      1000,...
    'bias',          1,...
    'uniquesv',      0,...
    'autosave',      600,...
    'filename',      'pca_f_autosave',...
    'minangle',      1e-8,...
    'algorithm',     'vb',...
    'niter_broadprior', 100,...
    'earlystop',     0,...
    'rmsstop',       [ 100 1e-4 1e-3 ],... % [] means no rms stop criteria
    'cfstop',        [],... % [] means no cost stop criteria
    'verbose',       1,...
    'xprobe',        [],...
    'rotate2pca',    1,...
    'display',       0 );

[ opts, errmsg, wrnmsg ] = argschk( opts, varargin{:} );
if ~isempty(errmsg), error( errmsg ), end
if ~isempty(wrnmsg), warning( wrnmsg ), end
Xprobe = opts.xprobe;

switch opts.algorithm
case 'ppca'
    use_prior = 0;
    use_postvar = 0;
case 'map'
    use_prior = 1;
    use_postvar = 0;
case 'vb'
    use_prior = 1;
    use_postvar = 1;
otherwise
    error( 'Wrong value of the argument ''algorithm''' )
end

[n1x,n2x] = size(X);
[ X, Xprobe, Ir, Ic, opts.init ] = rmempty( X, Xprobe,...
                                            opts.init, opts.verbose );

[n1,n2] = size(X);

if issparse(X)
    % X is a sparse matrix with only observed values
    M = spones(X);    Mprobe = spones(Xprobe);
else
    % Missing values are marked as NaNs
    M = ~isnan(X);    Mprobe = ~isnan(Xprobe);

    X(X==0) = eps;    Xprobe(Xprobe==0) = eps;
    X(isnan(X)) = 0;  Xprobe(isnan(Xprobe)) = 0;
    
end
Nobs_i = full(sum(M,2));
ndata = sum(Nobs_i);

nprobe = nnz(Mprobe);
if nprobe == 0
    Xprobe = [];
    opts.earlystop = 0;
end
    
[IX,JX,data] = find(X);
clear data

% Compute indices Isv: Sv{Isv(j)} gives Sv for j, j=1...n2
if opts.uniquesv
    [ nobscomb, obscombj, Isv ] = miscomb(M,opts.verbose);
else
    nobscomb = n2; Isv = []; obscombj = {};
end

[ A, S, Mu, V, Av, Sv, Muv ] = ...
    InitParms( opts.init, n1, n2, ncomp, nobscomb, Isv );

if use_prior
    Va = 1000*ones(1,ncomp); Vmu = 1000;
else
    Va = repmat(inf,1,ncomp); Vmu = inf;
end

% MAP estimation for Mu and A
if ~use_postvar, Muv = []; Av = {}; end
if ~opts.bias, Muv = []; Vmu = 0; end

if isempty(Mu)
    if opts.bias
        Mu = sum(X,2) ./ Nobs_i;
    else
        Mu = zeros(n1,1);
    end
end
[X,Xprobe] = SubtractMu( Mu, X, M, Xprobe, Mprobe, opts.bias );

[rms,errMx] = compute_rms( X, A, S, M, ndata );
prms = compute_rms( Xprobe, A, S, Mprobe, nprobe );

lc.rms = rms; lc.prms = prms; lc.time = 0; lc.cost = NaN;

dsph = DisplayInit( opts.display, lc );
PrintFirstStep( opts.verbose, rms, prms );
Aold = A;

% Parameters of the prior for variance parameters
hpVa = 0.001; hpVb = 0.001; hpV = 0.001;

time_start = clock;
time_autosave = time_start;
tic
%opts.niter_broadprior = 100;
for iter = 1:opts.maxiters

    % The prior is not updated at the begining of learning
    % to avoid killing sources 
    if use_prior && iter > opts.niter_broadprior
        % Update Va, Vmu
        
        if opts.bias
            Vmu = sum( Mu.^2 );
            if ~isempty(Muv)
                Vmu = Vmu + sum(Muv);
            end
            Vmu = (Vmu + 2*hpVa) / (n1 + 2*hpVb);
        end
        
        Va = sum( A.^2, 1 );
        if ~isempty(Av)
            for i = 1:n1
                Va = Va + diag(Av{i})';
            end
        end
        Va = (Va + 2*hpVa) / (n1 + 2*hpVb);
    end
    
    if opts.bias
        dMu = full( sum(errMx,2) ./ Nobs_i );
        if ~isempty(Muv)
            Muv = V ./ ( Nobs_i + V/Vmu );
        end
        th = 1 ./ ( 1 + V./Nobs_i/Vmu );  
        Mu_old = Mu;
        Mu = th.*( Mu + dMu );
        dMu = Mu - Mu_old;
        [X,Xprobe] = SubtractMu( dMu, X, M, Xprobe, Mprobe, 1 );
    end
    
    % Update S
    if isempty(Isv)
        for j = 1:n2
            %A_j = repmat(full(M(:,j)),1,ncomp) .* A;
            A_j = repmat(M(:,j),1,ncomp) .* A;
            Psi = A_j' * A_j + diag( repmat(V,1,ncomp) );
            if ~isempty(Av)
                for i = find(M(:,j))'
                    Psi = Psi + Av{i};
                end
            end
            invPsi = inv(Psi);
            S(:,j) = invPsi * A_j' * X(:,j);
            Sv{j} = V * invPsi;

            PrintProgress( opts.verbose, j, n2, 'Updating S:' )
        end
    else
        for k = 1:nobscomb
            j = obscombj{k}(1);
            %A_j = repmat(full(M(:,j)),1,ncomp) .* A;
            A_j = repmat(M(:,j),1,ncomp) .* A;
            Psi = A_j' * A_j + diag( repmat(V,1,ncomp) );
            if ~isempty(Av)
                for i = find(M(:,j))'
                    Psi = Psi + Av{i};
                end
            end
            invPsi = inv(Psi);
            Sv{k} = V * invPsi;
            tmp = invPsi * A_j';
            for j = obscombj{k}
                S(:,j) = tmp * X(:,j);
            end
            %S(:,obscombj{k}) = tmp * X(:,obscombj{k});
            PrintProgress( opts.verbose, k, nobscomb, 'Updating S:' )
        end
    end
    if opts.verbose == 2, fprintf('\r'), end

    if opts.rotate2pca
        [ dMu, A, Av, S, Sv ] = RotateToPCA( ...
            A, Av, S, Sv, Isv, obscombj, opts.bias );
        if opts.bias, 
            [X,Xprobe] = SubtractMu( dMu, X, M, Xprobe, Mprobe, 1 );
            Mu = Mu + dMu;
        end
    end
    
    % Update A
    if opts.verbose == 2
        fprintf('                                              \r')
    end
    for i = 1:n1
        %S_i = repmat(full(M(i,:)),ncomp,1) .* S;
        S_i = repmat(M(i,:),ncomp,1) .* S;
        Phi = S_i * S_i' + diag(V./Va);
        for j = find(M(i,:))
            if isempty(Isv)
                Phi = Phi + Sv{j};
            else 
                Phi = Phi + Sv{Isv(j)};
            end
        end
        invPhi = inv(Phi);
        A(i,:) = X(i,:) * S_i' * invPhi;
        if ~isempty(Av)
            Av{i} = V * invPhi;
        end
        
        PrintProgress( opts.verbose, i, n1, 'Updating A:' )
    end
    if opts.verbose == 2, fprintf('\r'), end
    
    [rms,errMx] = compute_rms( X, A, S, M, ndata );
    prms = compute_rms( Xprobe, A, S, Mprobe, nprobe );
    
    % Update V
    sXv = 0;
    if isempty(Isv)
        for r = 1:ndata
            sXv = sXv + A(IX(r),:) * Sv{JX(r)} * A(IX(r),:)';
            if ~isempty(Av)
                sXv = sXv + S(:,JX(r))' * Av{IX(r)} * S(:,JX(r)) + ...
                    sum( sum( Sv{JX(r)} .* Av{IX(r)} ) );
            end
        end
    else
        for r = 1:ndata
            sXv = sXv + A(IX(r),:) * Sv{Isv(JX(r))} * A(IX(r),:)';
            if ~isempty(Av)
                sXv = sXv + ...
                    S(:,JX(r))' * Av{IX(r)} * S(:,JX(r)) + ...
                    sum( sum( Sv{Isv(JX(r))} .* Av{IX(r)} ) );
            end
        end
    end

    if ~isempty(Muv)
        sXv = sXv + sum(Muv(IX));
    end
    
    sXv = sXv + (rms^2)*ndata;
    
    %V = rms^2 + V/ndata;
    V = ( sXv + 2*hpV ) / (ndata + 2*hpV);
    
    t = toc;
    lc.rms = [ lc.rms rms ]; lc.prms = [ lc.prms prms ];
    lc.time = [ lc.time t ];
    
    if ~isempty(opts.cfstop)
        cost = ...
            cf_full( X, A, S, Mu, V, Av, Sv, Isv, Muv, Va, Vmu,...
                      M, sXv, ndata );
        %cost = cf_full( X, A, S, Mu, V, Va, Vmu, M, sXv, ndata );
        %cost = cf_ppca( X, M, A.me, S.me, Mu, V, S.va, S.Iv );
        lc.cost = [ lc.cost cost ];
    end
    
    DisplayProgress( dsph, lc )
    angleA = subspace(A,Aold);
    PrintStep( opts.verbose, lc, angleA )

    convmsg = converg_check( opts, lc, angleA );
    if ~isempty(convmsg)
        if use_prior && iter <= opts.niter_broadprior
            % if the prior has never been updated: do nothing
        elseif opts.verbose, fprintf( '%s', convmsg ), end
        break
    end
    Aold = A;

    time = clock;
    if etime(time,time_autosave) > opts.autosave
        time_autosave = time;
        if opts.verbose == 2, fprintf('Saving ... '), end
        save( opts.filename,...
              'A', 'S', 'Mu', 'V', 'Av', 'Muv', 'Sv', 'Isv',...
              'Va', 'Vmu',...
              'lc', 'Ir', 'Ic', 'n1x', 'n2x', 'n1', 'n2' )
        if opts.verbose == 2, fprintf('done\n'), end
    end
end

% Finally rotate to the PCA solution
if ~opts.rotate2pca
    [ dMu, A, Av, S, Sv ] = RotateToPCA( ...
        A, Av, S, Sv, Isv, obscombj, opts.bias );
    if opts.bias, 
        Mu = Mu + dMu;
    end
end

if n1 < n1x
    [ A, Av ] = addmrows( A, Av, Ir, n1x, Va );
    [ Mu, Muv ] = addmrows( Mu, Muv, Ir, n1x, Vmu );
end
if n2 < n2x
    [ S, Sv, Isv ] = addmcols( S, Sv, Ic, n2x, Isv );
end

if nargout == 1
    A = struct( 'A', A );
    A.S = S;
    A.Mu = Mu;
    A.V = V;
    A.Va = Va;
    A.Vmu = Vmu;
    A.Av = Av;
    A.Sv = Sv;
    A.Isv = Isv;
    A.Muv = Muv;
    A.lc = lc;

else
    cv.A = Av;
    cv.S = Sv;
    cv.Isv = Isv;
    cv.Mu = Muv;

    hp.Va = Va;
    hp.Vmu = Vmu;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [X,Xprobe] = SubtractMu( Mu, X, M, Xprobe, Mprobe, update_bias )

n2 = size(X,2);

if ~update_bias
    return
end
   
if issparse(X)
    X = subtract_mu( X, Mu );
    if ~isempty(Xprobe)
        Xprobe = subtract_mu( Xprobe, Mu );
    end
else
    X = X - repmat(Mu,1,n2).*M;
    if ~isempty(Xprobe)
        Xprobe = Xprobe - repmat( Mu, 1, n2 ).*Mprobe;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find the PCA rotation: This has to be checked
function [ dMu, A, Av, S, Sv ] = ...
    RotateToPCA( A, Av, S, Sv, Isv, obscombj, update_bias );

n1 = size(A,1);
n2 = size(S,2);

if update_bias
    mS = mean(S,2);
    dMu = A*mS;
    S = S - repmat(mS,1,n2);
else
    dMu = 0;
end

covS = S*S';
if isempty(Isv)
    for j = 1:n2
        covS = covS + Sv{j};
    end
else
    nobscomb = length(obscombj);
    for j = 1:nobscomb
        covS = covS + ( length(obscombj{j})*Sv{j} );
    end
end
    
covS = covS / n2;
%covS = covS / (n2-n1);
[VS,D] = eig(covS);
RA = VS*sqrt(D);
A = A*RA;
covA = A'*A;
if ~isempty(Av)
    for i = 1:n1
        Av{i} = RA'*Av{i}*RA;
        covA = covA + Av{i};
    end
end
covA = covA / n1;
[VA,DA] = eig(covA);
[DA,I] = sort( -diag(DA) );
DA = -DA;
VA = VA(:,I);
A = A*VA;

if ~isempty(Av)
    for i = 1:n1
        Av{i} = VA'*Av{i}*VA;
    end
end
R = VA'*diag(1./sqrt(diag(D)))*VS';

S = R*S;
for j = 1:length(Sv)
    Sv{j} = R*Sv{j}*R';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ A, S, Mu, V, Av, Sv, Muv ] = ...
    InitParms( init, n1, n2, ncomp, nobscomb, Isv )

if ischar(init)
    if strcmpi( init, 'random' )
        init = struct([]);
    else
        % Load from a file
        init = load( init );
    end
end

if isstruct(init)
    if isfield( init, 'A' )
        A = init.A;
    else
        A = orth(randn(n1,ncomp));
    end
    if isfield( init, 'Av' ) && ~isempty(init.Av)
        if iscell(init.Av)
            Av = init.Av;
        else
            for i = 1:n1, Av{i} = diag(init.Av(i,:)); end
        end
    else
        Av = cell(1,n1);
        for i = 1:n1
            Av{i} = eye(ncomp);
        end
    end

    if isfield( init, 'Mu' )
        Mu = init.Mu;
    else
        Mu = [];
    end
    if isfield( init, 'Muv' ) && ~isempty(init.Muv)
        Muv = init.Muv;
    else
        Muv = ones(n1,1);
    end

    if isfield( init, 'V' )
        V = init.V;
    else
        V = 1;
    end
    
    if isfield( init, 'S' )
        S = init.S;
    else
        S = randn(ncomp,n2);
    end

    if isfield( init, 'Sv' ) && ~isempty(init.Sv)
       if nobscomb < n2
            [B,I] = unique(Isv,'first');
            if ~iscell(init.Sv)
                Sv = cell(1,nobscomb);
                for j = 1:nobscomb
                    Sv{j} = diag(init.Sv(:,Isv(I(j))));
                end
            
            elseif isfield( init, 'Isv' ) && ~isempty(init.Isv)
                Sv = { init.Sv{ init.Isv(I) } };
            else
                for j = 1:nobscomb
                    Sv{j} = init.Sv{Isv(I(j))};
                end
            end
        else
            if ~iscell(init.Sv)
                Sv = cell(1,n2);
                for j = 1:n2, Sv{j} = diag(init.Sv(:,j)); end
            
            elseif isfield( init, 'Isv' ) && ~isempty(init.Isv)
                Sv = { init.Sv{ init.Isv } };
            elseif length(init.Sv) == n2
                Sv = init.Sv;
            end
       end
    else
        Sv = cell(1,nobscomb);
        for j = 1:nobscomb, Sv{j} = eye(ncomp); end
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function PrintFirstStep( verbose, rms, prms )

if ~verbose
    return
end
fprintf( 'Step 0: rms = %.6f', rms )
if ~isnan(prms)
    fprintf( ' (%.6f)', prms )
end
fprintf( '\n' )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function PrintStep( verbose, lc, Aangle )

if ~verbose
    return
end

iter = length(lc.rms)-1;
steptime = lc.time(end)-lc.time(end-1);

fprintf( 'Step %i: ', iter )
if ~isnan(lc.cost(end))
    fprintf( 'cost = %.6f, ', lc.cost(end) );
end
fprintf( 'rms = %.6f', lc.rms(end) );
if ~isnan(lc.prms(end))
    fprintf( ' (%.6f)', lc.prms(end) );
end
fprintf( ', angle = %.2e', Aangle )
if steptime > 1
    fprintf( ' (%i sec)\n', round(steptime) )
else
    fprintf( ' (%.0e sec)\n', steptime )
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function PrintProgressBar( verbose, str )

if verbose == 2
    fprintf( [ str '\n' ] )
    %fprintf( '\n|                                                  |\r|' )
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function PrintProgress( verbose, i, n, str )

if verbose == 2
    fprintf( '\r%s %i/%i', str, i, n )
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dsph = DisplayInit( display, lc )

dsph.display = display;
if ~dsph.display
    return
end
dsph.fig = figure;
subplot(2,1,1)
dsph.rms = plot( 0:length(lc.rms)-1, lc.rms );
title( 'RMS training error' )
subplot(2,1,2)
title( 'RMS test error' )
dsph.prms = plot( 0:length(lc.prms)-1, lc.prms );

drawnow

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function DisplayProgress( dsph, lc )

if dsph.display
    set( dsph.rms, 'xdata', 0:length(lc.rms)-1,...
                   'ydata', lc.rms )
    set( dsph.prms, 'xdata', 0:length(lc.prms)-1,...
                   'ydata', lc.prms )
    drawnow
end
