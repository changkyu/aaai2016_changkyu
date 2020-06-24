%  PCA_PT - PCA with maximum a posteriori estimation
%
%  PCA with maximum a posteriori estimates for the parameters of the
%  model X(:,j) = Mu + A*S(:,j) + Noise. The noise is isotropic and
%  Gaussian with the variance V. Parameters Mu, A and S are assigned
%  Gaussian priors.
%
%  [ A, S, Mu, V, HP, LC ] = PCA_PT( X, N ) identifies the model for
%  the given data matrix X and number of principal components N.
%  Matrix X can be either sparse with only observed values included
%  (note that observed zeros should be replaced with eps) or a full
%  matrix with missing values replaced by NaNs.
%
%  The function returns the mean values of the model parameters in A,
%  S, Mu, and V. HP contains point estimates of the hyperparameters:
%  HP.Va is the prior variance for A(:,k) and HP.Vmu is the same for
%  Mu.
%  
%  LC is a structure with learning curves (rms training and probing
%  errors, cost function values, time).
%
%  Other optional parameter/value pairs with {default values}:
%   init       - Initialization type ({'random'} - random,
%                filename: load from file, structure: from given data)
%   gradient_type - Type of the gradient procedure used ('standard' for
%                gradient descent, {'alpha'} for modified Newton,
%                'minimize' for alternate update of A and S)
%   maxiters   - Maximum number of iterations {1000}
%   uniquesv   - Whether to use a speed-up trick for
%                gradient_type = 'minimize'
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
%  OUT = PCA_PT( X, N ) returns all the outputs in a single
%  structure variable OUT. Learning can be continued as follows:
%    out = pca_pt( X, N );
%    [ A, S, Mu, V, HP, LC ] = pca_pt( X, N, 'init', out );

%  This software is provided "as is", without warranty of any kind.
%  Alexander Ilin, Tapani Raiko

function [ A, S, Mu, V, hp, lc ] = pca_pt( X, ncomp, varargin )

% Default parameter values
defopts = struct( ...
    'gradient_type', 'alpha',...
    'init',          'random',...
    'minangle',      1e-8,...
    'maxiters',      inf,...
    'uniquesv',      1,...
    'bias',          1,...
    'autosave',      3600,...
    'earlystop',     0,...
    'rmsstop',       [],... % [] means no rms stop criteria
    'cfstop',        [ 100 1e-4 1e-3 ],... % [] means no cost stop criteria
    'verbose',       1,...
    'xprobe',        [],...
    'filename',      'pca_p_autosave',...
    'rotate2pca',    0,...
    'display',       0 );
numCPU = 1;

[ opts, errmsg, wrnmsg ] = argschk( defopts, varargin{:} );
if ~isempty(errmsg), error( errmsg ), end
if ~isempty(wrnmsg), warning( wrnmsg ), end
Xprobe = opts.xprobe;

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
end

if ~strcmpi( opts.gradient_type, 'minimize' )
    opts.uniquesv = 0;
end

% Compute the number of missing values combinations
if opts.uniquesv
    [ nobscomb, obscombj ] = miscomb(M,opts.verbose);
else
    nobscomb = n2; obscombj = {};
end

[ AA, SS, Mu, V ] = InitParms( opts.init, n1, n2, ncomp );

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

[rms,errMx] = compute_rms( X, AA{cur}, SS{cur}, M, ndata );
prms = compute_rms( Xprobe, AA{cur}, SS{cur}, Mprobe, nprobe );

lc.rms = rms; lc.prms = prms; lc.time = 0; lc.cost = NaN;

% Parameters of the prior for variance parameters
hpVa = 0.001; hpVb = 0.001; hpV = 0.001;
Va = ones(1,ncomp);
Vmu = 1;
%Va(:) = 1000;
%Vmu = 1000;
if ~opts.bias
    Vmu = 0;
end

time_start = clock;
time_autosave = time_start;
sd_iter = 0;
tic
for iter = 1:opts.maxiters

    % Update Va, Vmu
    %Va = mean( AA{cur}.^2, 1 );
    % Regularized: Inverse-gamma prior
    if iter > 10
        Va = sum( AA{cur}.^2 + 2*hpVa, 1 )/(n1 + 2*hpVb);
        if opts.bias
            Vmu = sum( Mu.^2 + 2*hpVa, 1 )/(n1 + 2*hpVb);  %mean( Mu.^2 );
        end
    end
    if opts.bias && iter > 1
        dMu = full( sum(errMx,2) ./ Nobs_i );
        th = 1 ./ ( 1 + V./Nobs_i/Vmu );   %th = Muv .* Nobs_i / V;
        Mu_old = Mu;
        Mu = th.*( Mu + dMu );
        dMu = Mu - Mu_old;
        [X,Xprobe] = SubtractMu( dMu, X, M, Xprobe, Mprobe, 1 );
    end

    % Update A, S
    if strcmp(opts.gradient_type,'minimize')
        
        if iter == 1 && opts.verbose
            fprintf( 'Step 0: rms=%.6f\n', rms )
        end

        % Update S
        if nobscomb == n2
            for j = 1:n2
                A_j = repmat(full(M(:,j)),1,ncomp) .* AA{cur};
                Psi = A_j' * A_j + diag( repmat(V,1,ncomp) );
                invPsi = inv(Psi);
                SS{new}(:,j) = invPsi * A_j' * X(:,j);
                
                PrintProgress( opts.verbose, j, n2, 'Updating S:' )
            end
            
        else
            for k = 1:nobscomb
                j = obscombj{k}(1);
                A_j = repmat(full(M(:,j)),1,ncomp) .* AA{cur};
                Psi = A_j' * A_j + diag( repmat(V,1,ncomp) );
                invPsi = inv(Psi);
                tmp = invPsi * A_j';
                for j = obscombj{k}
                    SS{new}(:,j) = tmp * X(:,j);
                end
                PrintProgress( opts.verbose, k, nobscomb, 'Updating S:' )
            end
        end
        if opts.verbose, fprintf('\r'), end
        
        if opts.rotate2pca
            [ dMu, AA{cur}, SS{new} ] = ...
                RotateToPCA( AA{cur}, SS{new}, opts.bias );
            if opts.bias, 
                [X,Xprobe] = SubtractMu( ...
                    dMu, X, M, Xprobe, Mprobe, 1 );
                Mu = Mu + dMu;
            end
        end
        
        % Update A
        if opts.verbose
            fprintf('                                              \r')
        end
        for i = 1:n1
            S_i = repmat(full(M(i,:)),ncomp,1) .* SS{new};
            Phi = S_i * S_i' + diag(V./Va);
            invPhi = inv(Phi);
            AA{new}(i,:) = X(i,:) * S_i' * invPhi;
        
            PrintProgress( opts.verbose, i, n1, 'Updating A:' )
        end
        if opts.verbose, fprintf('\r'), end

        if isempty(opts.cfstop)
            [rms,errMx] = compute_rms( X, AA{new}, SS{new}, M, ndata );
            cost = [];
        else
            [cost,errMx,rms] = cf_pt( ...
                X, AA{new}, SS{new}, Mu, V, Va, Vmu, M, ndata, numCPU );
            %[errMx,rms,cost] = cf_map( X, AA{new}, SS{new},...
            %                               Mu, Va, V, Vmu, M,...
            %                               ndata, numCPU );

        end
        
        prms = compute_rms( Xprobe, AA{new}, SS{new}, Mprobe, nprobe );
        
    else
        
        [cost_old,errMx,rms] = cf_pt( ...
            X, AA{cur}, SS{cur}, Mu, V, Va, Vmu, M, ndata, numCPU );
        %[errMx,rms,cost_old] = cf_map( X, AA{cur}, SS{cur},...
        %                               Mu, Va, V, Vmu, M,...
        %                               ndata, numCPU );
        
        if iter == 1
            lc.cost = cost_old;
            if opts.verbose
                fprintf( 'Step 0: cost=%.6e, rms=%.6f\n', cost_old, rms )
            end
        end
        
        switch opts.gradient_type
        case 'normal'
            AA{new} = AA{cur} + ...
                      full( SS{cur} * errMx' )'*(lrateA/V) ...
                      - AA{cur}.* repmat(lrateA./Va,n1,1);
            
            SS{new} = SS{cur} + full( ((lrateS/V) * AA{cur}') * errMx ) ...
                     - SS{cur} * lrateS;
            
        case 'alpha'
            dgHsA = ( repmat(1./Va,n1,1) + full( SS{cur}.^2 * M' )'/V );
            dgHsS = ( 1 + full( ((AA{cur}.^2)/V)' * M ) );
            
            AA{new} = AA{cur} + dgHsA.^(-alpha) .* ...
                     ( full( SS{cur} * errMx' )'*(lrateA/V) ...
                       - AA{cur} .* repmat(lrateA./Va,n1,1) );
        
            SS{new} = SS{cur} + dgHsS.^(-alpha) .* ...
                     (full( ((lrateS/V) * AA{cur}') * errMx ) ...
                      - SS{cur} * lrateS);
        end
    
        [cost,errMx,rms] = cf_pt( ...
            X, AA{new}, SS{new}, Mu, V, Va, Vmu, M, ndata, numCPU );
        %[errMx,rms,cost] = cf_map( X, AA{new}, SS{new},...
        %                           Mu, Va, V, Vmu, M,...
        %                           ndata, numCPU );
        
        if cost < cost_old
            lrate = lrate * 1.2;
            lrateA = lrateA * 1.2;
            lrateS = lrateS * 1.2;
        else
            if opts.verbose
                fprintf('Slowing down (%e) ', cost_old-cost );
            end
            sd_iter = 0;
            while cost > cost_old && sd_iter < 10
                if opts.verbose, fprintf('.'); end
                lrate = lrate * 0.5;
                lrateA = lrateA * 0.5;
                lrateS = lrateS * 0.5;
                AA{new} = (AA{new}+AA{cur})/2;
                SS{new} = (SS{new}+SS{cur})/2;
                
                [cost,errMx,rms] = cf_pt( ...
                    X, AA{new}, SS{new}, Mu, V, Va, Vmu,...
                    M, ndata, numCPU );
                %[errMx,rms,cost] = cf_map( X, AA{new}, SS{new},...
                %                           Mu, Va, V, Vmu, M,...
                %                           ndata, numCPU );
                sd_iter = sd_iter + 1;
            end
            if opts.verbose, fprintf(' ok\n'); end
        end
        
        if opts.rotate2pca
            [ dMu, AA{new}, SS{new} ] = ...
                RotateToPCA( AA{new}, SS{new}, opts.bias );
            if opts.bias, 
                [X,Xprobe] = SubtractMu( ...
                    dMu, X, M, Xprobe, Mprobe, opts.bias );
                Mu = Mu + dMu;
            end

            Va = sum( AA{cur}.^2 + 2*hpVa, 1 )/(n1 + 2*hpVb);
            if opts.bias
                Vmu = sum( Mu.^2 + 2*hpVa, 1 )/(n1 + 2*hpVb);  %mean( Mu.^2 );
            end
        end

        prms = compute_rms( Xprobe, AA{new}, SS{new}, Mprobe, nprobe );
        
    end

    %V = rms^2;
    V = ( rms^2*ndata + 2*hpV ) / (ndata + 2*hpV);

    t = toc;
    lc.rms = [ lc.rms rms ]; lc.prms = [ lc.prms prms ];
    lc.cost = [ lc.cost cost ]; lc.time = [ lc.time t ];

    %DisplayProgress( dsph, lc )
    angleA = subspace(AA{new},AA{cur});
    PrintStep( opts.verbose, lc, angleA )

    convmsg = converg_check( opts, lc, angleA, sd_iter );
    if ~isempty(convmsg)
        if opts.verbose, fprintf( '%s', convmsg ), end
        break
    end
    
    time = clock;
    if etime(time,time_autosave) > opts.autosave
        if opts.verbose==2, fprintf( 'Saving ... ' ), end
        time_autosave = time;
        A = AA{new};
        S = SS{new};
        save( opts.filename,...
              'A', 'S', 'Mu', 'V',...
              'lc', 'Ir', 'Ic', 'n1x', 'n2x', 'n1', 'n2' )
        clear A S
        if opts.verbose==2, fprintf( 'done\n' ), end
    end

    if cur == 1,  cur = 2; new = 1;
    else          cur = 1; new = 2;
    end

end

A = AA{cur};
S = SS{cur};
clear AA SS

% Sorting according to explained variance should be done by RotateToPCA

if n1 < n1x
    A = addmrows( A, [], Ir, n1x, Va );
    Mu = addmrows( Mu, [], Ir, n1x, Vmu );
end
if n2 < n2x
    S = addmcols( S, [], Ic, n2x );
end

if nargout == 1
    A = struct( 'A', A );
    A.S = S;
    A.Mu = Mu;
    A.V = V;
    A.Va = Va;
    A.Vmu = Vmu;
    A.lc = lc;

else
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
function [ dMu, A, S ] = RotateToPCA( A, S, update_bias )

n1 = size(A,1);
n2 = size(S,2);

if update_bias
    mS = mean(S,2);
    dMu = A*mS;
    S = S - repmat(mS,1,n2);
else
    dMu = 0;
end

covS = S*S' / n2;
[VS,D] = eig(covS);
RA = VS*sqrt(D);
A = A*RA;

covA = A'*A / n1;
[VA,DA] = eig(covA);
[DA,I] = sort( -diag(DA) );
DA = -DA;
VA = VA(:,I);
A = A*VA;

R = VA'*diag(1./sqrt(diag(D)))*VS';

S = R*S;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ AA, SS, Mu, V ] = InitParms( init, n1, n2, ncomp )

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

    if isfield( init, 'Mu' )
        Mu = init.Mu;
    else
        Mu = [];
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

end

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
function PrintProgress( verbose, i, n, str )

if verbose == 2
    fprintf( '\r%s %i/%i', str, i, n )
end

