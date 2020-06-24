%  PCA_DIAG - PCA with fully factorial Gaussian posterior
%
%  PCA with fully factorial Gaussian pdf (diagonal covariance matrix)
%  for approximating the posterior distributions of A and S in the
%  model X(:,j) = Mu + A*S(:,j) + Noise. The noise is isotropic and
%  Gaussian with the variance V.
%
%  [ A, S, Mu, V, CV, HP, LC ] = PCA_DIAG( X, N ) identifies the model
%  for the given data matrix X and number of principal components N.
%  Matrix X can be either sparse with only observed values included
%  (note that observed zeros should be replaced with eps) or a full
%  matrix with missing values replaced by NaNs.
%
%  The function returns the mean values of the model parameters in A,
%  S, Mu, and V. The posterior variances are stored in CV such that
%  CV.A(i,k) is the posterior variance of A(i,k) and CV.S(k,j) is the
%  same for S(k,j). HP contains point estimates of the
%  hyperparameters: HP.Va is the prior variance for A(:,k) and HP.Vmu
%  is the same for Mu.
%  
%  LC is a structure with learning curves (rms training and probing
%  errors, cost function values, time).
%
%  PCA_DIAG( X, N, 'algorithm', name ) specifies the algorithm:
%   'ppca': Probabilistic PCA (no prior and point estimates for A, Mu)
%   'map':  Prior on A(:,k) and Mu, point estimates for A(i,:) and Mu
%   'vb':   Variational Bayesian PCA (prior on A and Mu, Gaussian
%           posterior approximation for A(i,:) and Mu) {default}
%  
%  Other optional parameter/value pairs with {default values}:
%   init       - Initialization type ({'random'} - random,
%                filename: load from file, structure: from given data)
%   gradient_type - Type of the gradient procedure used ('standard' for
%                gradient descent, {'alpha'} for modified Newton)
%   maxiters   - Maximum number of iterations {1000}
%   numcpu     - Number of CPUs for parallel computing (only for
%                sparse data sets)
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
%   autosave   - Auto-save after each {3600} seconds
%   savebest     Whether to save the model with the best probing error
%   filename   - Name of the file for auto-save {'pca_d_autosave'}
%
%  OUT = PCA_DIAG( X, N ) returns all the outputs in a single
%  structure variable OUT. Learning can be continued as follows:
%    out = pca_diag( X, N );
%    [ A, S, Mu, V, CV, HP, LC ] = pca_diag( X, N, 'init', out );

%  TODO: Order the output components based on the variance.

%  This software is provided "as is", without warranty of any kind.
%  Alexander Ilin, Tapani Raiko

function [ A, S, Mu, V, cv, hp, lc ] = pca_diag( X, ncomp, varargin )

% Default parameter values
defopts = struct( ...
    'gradient_type', 'alpha',...
    'init',          'random',...
    'minangle',      1e-8,...
    'maxiters',      inf,...
    'algorithm',     'vb',...
    'earlystop',     0,...
    'niter_broadprior', 100,...
    'rmsstop',       [],... % [] means no cost stop criteria
    'cfstop',        [ 100 1e-4 1e-3 ],... % [] means no cost stop criteria
    'verbose',       1,...
    'bias',          1,...
    'autosave',      3600,...
    'savebest',      0,...   % Save A and S with the best prms
    'xprobe',        [],...
    'filename',      'pca_d_autosave',...
    'numcpu',        1,...
    'display',       0 );

[ opts, errmsg, wrnmsg ] = argschk( defopts, varargin{:} );
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
Nobs_i = sum(M,2);
ndata = sum(Nobs_i);

nprobe = nnz(Mprobe);
if nprobe == 0
    Xprobe = [];
    opts.earlystop = 0;
    if opts.savebest
        fprintf( 'No test data given for saving the best prms solution.\n' );
        opts.savebest = 0;
    end
end

[ AA, SS, Mu, V, Av, Sv, Muv ] = InitParms( opts.init, n1, n2, ncomp );

if use_prior
    Va = 1000*ones(1,ncomp); Vmu = 1000;
else
    Va = repmat(inf,1,ncomp);  Vmu = inf;
end

% MAP estimation for Mu and A
if ~use_postvar, Muv = []; Av = []; sumAv = 0; end
if ~opts.bias, Muv = []; Vmu = 0; end

if isempty(Mu)
    if opts.bias
        Mu = sum(X,2) ./ Nobs_i;
    else
        Mu = zeros(n1,1);
    end
end
[X,Xprobe] = SubtractMu( Mu, X, M, Xprobe, Mprobe, opts.bias );

cur = 1; new = 2;
lrate = 1; lrateA = 1; lrateS = 1;
alpha = 2/3;

prms = compute_rms( Xprobe, AA{cur}, SS{cur}, Mprobe, nprobe );

% Parameters of the prior for variance parameters
hpVa = 0.001; hpVb = 0.001; hpV = 0.001;

time_start = clock;
time_autosave = time_start;
sd_iter = 0;
best.prms = inf;

sumSv = full( Sv * M' )'/V;

tic
lc.time = 0;

for iter = 1:opts.maxiters

    % Update Va, Vmu
    if use_prior && iter > opts.niter_broadprior
        if isempty(Muv)
            Vmu = sum( Mu.^2 );
        else
            Vmu = sum( Mu.^2 + Muv );
        end
        Vmu = (Vmu + 2*hpVa) / (n1 + 2*hpVb);

        if isempty(Av)
            Va = sum( AA{cur}.^2, 1 );
        else
            Va = sum( AA{cur}.^2 + Av, 1 );
        end
        Va = (Va + 2*hpVa) / (n1 + 2*hpVb);
    end

    if opts.bias && iter > 1
        dMu = full( sum(errMx,2) ./ Nobs_i );
        if ~isempty(Muv)
            Muv = V ./ ( Nobs_i + V/Vmu );
        end
        th = 1 ./ ( 1 + V./Nobs_i/Vmu );   %th = Muv .* Nobs_i / V;
        Mu_old = Mu;
        Mu = th.*( Mu + dMu );
        dMu = Mu - Mu_old;
        [X,Xprobe] = SubtractMu( dMu, X, M, Xprobe, Mprobe, 1 );
    end

    if ~isempty(Av) || ~strcmpi(opts.gradient_type,'standard')
        Av = ( repmat(1./Va,n1,1) + sumSv + ...
               full( SS{cur}.^2 * M' )'/V ).^(-1);
        sumAv = full( (Av/V)' * M );
    end        

    Sv = ( 1 + sumAv + full( ((AA{cur}.^2)/V)' * M ) ).^(-1);
    sumSv = full( Sv * M' )'/V;

    if ~use_postvar
        [cost_old,errMx,rms] = ...
            cf_diag( X, AA{cur}, SS{cur}, Mu, V, [], Sv, [], Va, Vmu, M,...
                     Nobs_i, ndata, opts.numcpu );
        
    else
        [cost_old,errMx,rms] = ...
            cf_diag( X, AA{cur}, SS{cur}, Mu, V, Av, Sv, Muv, Va, Vmu, M,...
                     Nobs_i, ndata, opts.numcpu );
    end
    
    if iter == 1
        lc.rms = rms; lc.prms = prms; lc.cost = cost_old;
        if opts.verbose
            fprintf( 'Step 0: cost=%.6e, rms=%.6f\n', cost_old, rms )
        end
    end

    % Update A, S
    switch opts.gradient_type
    case 'standard'
        AA{new} = AA{cur} + ...
                  ( full( SS{cur} * errMx' )'/V - ...
                    AA{cur} .* sumSv - ...
                    AA{cur}.* repmat(1./Va,n1,1) ) ...
                  * lrateA;
        
        SS{new} = SS{cur} + ...
                  full( ((lrateS/V) * AA{cur}') * errMx ) -...
                  SS{cur} .* sumAv * lrateS - ...
                  SS{cur} * lrateS;
        
    case 'alpha'
        AA{new} = AA{cur} + Av.^alpha .* ...
                  ( full( SS{cur} * errMx' )'/V - ...
                    AA{cur} .* sumSv - ...
                    AA{cur} .* repmat(1./Va,n1,1) ) ...
                  * lrateA;
        
        SS{new} = SS{cur} + Sv.^alpha .* ...
                  (full( ((lrateS/V) * AA{cur}') * errMx ) - ...
                   SS{cur} .* sumAv * lrateS - ...
                   SS{cur} * lrateS);
        if ~use_postvar
            Av = [];
        end
    end
    
    [cost,errMx,rms,vN] = cf_diag( ...
        X, AA{new}, SS{new}, Mu, V,  Av, Sv, Muv, Va, Vmu, M,...
        Nobs_i, ndata, opts.numcpu );
    
    if cost < cost_old
        lrate = lrate * 1.1;
        lrateA = lrateA * 1.1;
        lrateS = lrateS * 1.1;
    else
        if opts.verbose
            fprintf('Slowing down (%e) ', cost_old-cost );
        end
        sd_iter = 0;
        while cost > cost_old && sd_iter < 40
            if opts.verbose, fprintf('.'); end
            lrate = lrate * 0.5;
            lrateA = lrateA * 0.5;
            lrateS = lrateS * 0.5;
            AA{new} = (AA{new}+AA{cur})/2;
            SS{new} = (SS{new}+SS{cur})/2;
            
            [cost,errMx,rms,vN] = cf_diag( ...
                X, AA{new}, SS{new}, Mu, V, Av, Sv, Muv, Va, Vmu, M,...
                Nobs_i, ndata, opts.numcpu );
            sd_iter = sd_iter + 1;
        end
        if opts.verbose, fprintf(' ok\n'); end
    end
    
    prms = compute_rms( Xprobe, AA{new}, SS{new}, Mprobe, nprobe );

    % Update V
    %V = vN/ndata;
    V = ( vN + 2*hpV ) / (ndata + 2*hpV);
    
    t = toc;
    lc.rms = [ lc.rms rms ]; lc.prms = [ lc.prms prms ];
    lc.cost = [ lc.cost cost ]; lc.time = [ lc.time t ];

    if opts.savebest && lc.prms(end) < best.prms
        best.A = AA{new}; best.S = SS{new}; best.Mu = Mu;
        best.V = V; best.rms = rms; best.prms = prms;
    end
    
    %DisplayProgress( dsph, lc )
    angleA = subspace(AA{new},AA{cur});
    PrintStep( opts.verbose, lc, angleA )

    convmsg = converg_check( opts, lc, angleA, sd_iter );
    if ~isempty(convmsg)
        if use_prior && iter <= opts.niter_broadprior
            % if the prior has never been updated: do nothing
        elseif opts.verbose, fprintf( '%s', convmsg ), end
        break
    end

    time = clock;
    if etime(time,time_autosave) > opts.autosave
        if opts.verbose==2, fprintf( 'Saving ... ' ), end
        time_autosave = time;
        A = AA{new}; S = SS{new};
        save( opts.filename,...
              'A', 'S', 'Mu', 'Av', 'Muv', 'V', 'Va', 'Vmu', 'Sv',...
              'lc', 'Ir', 'Ic', 'n1x', 'n2x', 'n1', 'n2' )
        clear A S
        if opts.savebest
            save( opts.filename, 'best', '-append' )
        end
        if opts.verbose==2, fprintf( 'done\n' ), end
    end

    if cur == 1,  cur = 2; new = 1;
    else          cur = 1; new = 2;
    end

end

A = AA{cur};
S = SS{cur};

% Sort components according to explained variance
if ~isempty(Av)
    normA = sum( A.^2 + Av, 1 );
else
    normA = sum( A.^2, 1 );
end
[tmp,I] = sort(-normA);
A = A(:,I);
if ~isempty(Av)
    Av = Av(:,I);
    Va = Va(I);
end
S = S(I,:);
Sv = Sv(I,:);

if n1 < n1x
    [ A, Av ] = addmrows( A, Av, Ir, n1x, Va );
    [ Mu, Muv ] = addmrows( Mu, Muv, Ir, n1x, Vmu );
end
if n2 < n2x
    [ S, Sv ] = addmcols( S, Sv, Ic, n2x );
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
    A.Muv = Muv;
    A.lc = lc;
else
    cv.A = Av;
    cv.S = Sv;
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
function [ AA, SS, Mu, V, Av, Sv, Muv ] = InitParms( init, n1, n2, ncomp )

AA = cell(2,1);
SS = cell(2,1);
cur = 1; new = 2;

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
        AA{cur} = init.A;
    else
        AA{cur} = orth(randn(n1,ncomp));
    end
    if isfield( init, 'Av' )
        Av = init.Av;
    else
        Av = ones(n1,ncomp);
    end

    if isfield( init, 'Mu' )
        Mu = init.Mu;
    else
        Mu = [];
    end
    if isfield( init, 'Muv' )
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
        SS{cur} = init.S;
    else
        SS{cur} = randn(ncomp,n2);
    end

    % TODO: if init.Isv exists
    if isfield( init, 'Sv' )
        if iscell(init.Sv)
            Sv = zeros(ncomp,n2);
            for j = 1:n2
                Sv(:,j) = diag(init.Sv{j});
            end
        else
            Sv = init.Sv;
        end
    else
        Sv = ones(ncomp,n2);
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function PrintStep( verbose, lc, Aangle )

if ~verbose
    return
end

iter = length(lc.rms)-1;
steptime = lc.time(end)-lc.time(end-1);

if ~isnan(lc.prms(end))
    fprintf( ...
        'Step %i: cost=%.6e, rms=%.6f (%.6f), angle=%.2e (%d sec)\n', ...
        iter, lc.cost(end), lc.rms(end), lc.prms(end),...
        Aangle, round(steptime) )
else    
    fprintf( ...
        'Step %i: cost=%.6e, rms=%.6f, angle=%.2e (%d sec)\n', ...
             iter, lc.cost(end), lc.rms(end), Aangle, round(steptime) )
end


