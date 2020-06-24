%  LSPCA - Least-squares approach to PCA with missing values
%
%  [ A, S, Mu, LC ] = LSPCA( X, N ) finds the PCA decomposition X = AS
%  + [Mu .. Mu] for the given data matrix X and number of principal
%  components N using the unregularized least-squares approach. Matrix
%  X can be either sparse with only observed values included (note
%  that observed zeros should be replaced with eps) or a full matrix
%  with missing values replaced by NaNs.

%  LC is a structure with learning curves (rms training and probing
%  errors, time).
%
%  Optional parameter/value pairs with {default values}:
%   init       - Initialization type ({'random'} - random,
%                structure: from given data)
%   gradient_type - Type of the gradient procedure used ('standard' for
%                gradient descent, {'alpha'} for modified Newton,
%                'minimize' for alternate update of A and S)
%   maxiters   - Maximum number of iterations {1000}
%   rotate2pca - Whether to perform rotation of A and S to the PCA
%                basis after each iteration
%                convergence {1}
%   uniquesv   - Whether to use a speed-up trick for
%                gradient_type = 'minimize'
%   minangle   - Termination by minimum angle between subspaces
%                defined by A {1e-8}
%   rmsstop    - Termination by rms training error ([] - no rms stop)
%                {[ 100 1e-4 1e-3 ]}: 100 iterations, 1e-4 absolute
%                tolerance, 1e-3 relative tolerance,
%   xprobe     - Validation data set (of the same dimensions as X)
%   earlystop  - Whether to use early stopping based on probing error
%   verbose    - Progress information display level (0,{1},2)
%   display    - Plot progress {0}
%   autosave   - Auto-save after each {3600} seconds
%   filename   - Name of the file for auto-save {'lspca_autosave'}
%
%  OUT = LSPCA( X, N ) returns all the outputs in a single
%  structure variable OUT.

%  This software is provided "as is", without warranty of any kind.
%  Alexander Ilin, Tapani Raiko

function [ A, S, Mu, lc ] = lspca( X, ncomp, varargin )

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
    'verbose',       1,...
    'xprobe',        [],...
    'filename',      'lspca_autosave',...
    'rotate2pca',    1,...
    'display',       'on' );

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

% Compute the number of missing values combinations
if opts.uniquesv
    [ nobscomb, obscombj ] = miscomb(M,opts.verbose);
else
    nobscomb = n2; obscombj = {};
end

[ AA, SS, Mu ] = InitParms( opts.init, n1, n2, ncomp );
%keyboard

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
if opts.verbose
    fprintf( 'Step 0: rms=%.6f\n', rms )
end

lc.rms = rms; lc.prms = prms; lc.time = 0;

% Parameters of the prior for variance parameters
hpVa = 0.001; hpVb = 0.001; hpV = 0.001;
Va = ones(1,ncomp);
Vmu = 1;

time_start = clock;
time_autosave = time_start;
sd_iter = 0;
tic
for iter = 1:opts.maxiters

    if opts.bias && iter > 1
        dMu = full( sum(errMx,2) ./ Nobs_i );
        Mu = Mu + dMu;
        [X,Xprobe] = SubtractMu( dMu, X, M, Xprobe, Mprobe, 1 );
    end

    % Update A, S
    if strcmp(opts.gradient_type,'minimize')
        
        % Update S
        if nobscomb == n2
            for j = 1:n2
                A_j = repmat(full(M(:,j)),1,ncomp) .* AA{cur};
                %SS{new}(:,j) = inv(A_j'*A_j) * A_j' * X(:,j);
                SS{new}(:,j) = pinv(A_j) * X(:,j);
                PrintProgress( opts.verbose, j, n2, 'Updating S:' )
            end
            
        else
            for k = 1:nobscomb
                j = obscombj{k}(1);
                A_j = repmat(full(M(:,j)),1,ncomp) .* AA{cur};
                tmp = pinv(A_j);
                for j = obscombj{k}
                    SS{new}(:,j) = tmp * X(:,j);
                end
                PrintProgress( opts.verbose, k, nobscomb, 'Updating S:' )
            end
        end
        if opts.verbose, fprintf('\r'), end
        
        if opts.rotate2pca
            [ dMu, AA{new}, SS{new} ] = ...
                RotateToPCA( AA{cur}, SS{new}, opts.bias );
            [X,Xprobe] = SubtractMu( dMu, X, M, Xprobe, Mprobe, opts.bias );
            if opts.bias, Mu = Mu + dMu; end
        end
        
        % Update A
        if opts.verbose
            fprintf('                                              \r')
        end
        for i = 1:n1
            S_i = repmat(full(M(i,:)),ncomp,1) .* SS{new};
            %AA{new}(i,:) = X(i,:) * S_i' * inv(S_i*S_i');
            AA{new}(i,:) = X(i,:) * pinv(S_i);
        
            PrintProgress( opts.verbose, i, n1, 'Updating A:' )
        end
        if opts.verbose, fprintf('\r'), end

        [rms,errMx] = compute_rms( X, AA{new}, SS{new}, M, ndata );
        prms = compute_rms( Xprobe, AA{new}, SS{new}, Mprobe, nprobe );
    
    else
        
        rms_old = rms;
        switch opts.gradient_type
        case 'normal'
            AA{new} = AA{cur} + full( SS{cur} * errMx' )' * lrateA;
            SS{new} = SS{cur} + full( (lrateS * AA{cur}') * errMx );

        case 'alpha'
            dgHsA = full( SS{cur}.^2 * M' )';
            dgHsS = full( (AA{cur}.^2)' * M );
            
            AA{new} = AA{cur} + dgHsA.^(-alpha) .* ...
                     full( SS{cur} * errMx' )' * lrateA;
        
            SS{new} = SS{cur} + dgHsS.^(-alpha) .* ...
                      full( (lrateS * AA{cur}') * errMx );

        end
        
        [rms,errMx] = compute_rms( X, AA{new}, SS{new}, M, ndata );
        
        if rms < rms_old
            lrate = lrate * 1.2;
            lrateA = lrateA * 1.2;
            lrateS = lrateS * 1.2;
        else
            if opts.verbose
                fprintf('Slowing down (%e) ', rms_old-rms );
            end
            sd_iter = 0;
            while rms > rms_old && sd_iter < 10
                if opts.verbose, fprintf('.'); end
                lrate = lrate * 0.5;
                lrateA = lrateA * 0.5;
                lrateS = lrateS * 0.5;
                AA{new} = (AA{new}+AA{cur})/2;
                SS{new} = (SS{new}+SS{cur})/2;
                
                [rms,errMx] = compute_rms( X, AA{new}, SS{new}, M, ndata );
                sd_iter = sd_iter + 1;
            end
            if opts.verbose, fprintf(' ok\n'); end
        end

        if opts.rotate2pca
            [ dMu, AA{new}, SS{new} ] = ...
                RotateToPCA( AA{new}, SS{new}, opts.bias );
            [X,Xprobe] = SubtractMu( dMu, X, M, Xprobe, Mprobe, opts.bias );
            if opts.bias, Mu = Mu + dMu; end
        end

        prms = compute_rms( Xprobe, AA{new}, SS{new}, Mprobe, nprobe );
        
    end

    t = toc;
    lc.rms = [ lc.rms rms ]; lc.prms = [ lc.prms prms ];
    lc.time = [ lc.time t ];

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
        % fprintf( 'Saving ... ' )
        time_autosave = time;
        A = AA{new};
        S = SS{new};
        save( opts.filename,...
              'A', 'S', 'Mu',...
              'lc', 'Ir', 'Ic', 'n1x', 'n2x', 'n1', 'n2' )
        clear A S
        % fprintf( 'done\n' )
    end

    if cur == 1,  cur = 2; new = 1;
    else          cur = 1; new = 2;
    end

end

A = AA{cur};
S = SS{cur};
clear AA SS

% Finally rotate to the PCA solution
if ~opts.rotate2pca
    [ dMu, A, S ] = RotateToPCA( A, S, opts.bias );
    if opts.bias, 
        Mu = Mu + dMu;
    end
end

if n1 < n1x
    A = addmrows( A, [], Ir, n1x, NaN );
    Mu = addmrows( Mu, [], Ir, n1x, NaN );
end
if n2 < n2x
    S_new = repmat(NaN, ncomp, n2x);
    S_new(:,Ic) = S;
    S = S_new;
    clear S_new
end

if nargout == 1
    A = struct( 'A', A );
    A.S = S;
    A.Mu = Mu;
    A.lc = lc;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [X,Xprobe] = SubtractMu( Mu, X, M, Xprobe, Mprobe, update_bias )

n2 = size(X,2);

if update_bias
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
function [ AA, SS, Mu ] = InitParms( init, n1, n2, ncomp )

AA = cell(2,1);
SS = cell(2,1);
cur = 1; new = 2;

if ischar(init) && strcmpi( init, 'random' )
    init = struct([]);
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

if ~isnan(lc.prms(end))
    fprintf( ...
        'Step %i: rms=%.6f (%.6f), angle=%.2e (%d sec)\n', ...
        iter, lc.rms(end), lc.prms(end), Aangle, round(steptime) )
else    
    fprintf( ...
        'Step %i: rms=%.6f, angle=%.2e (%d sec)\n', ...
             iter, lc.rms(end), Aangle, round(steptime) )
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function PrintProgress( verbose, i, n, str )

if verbose == 2
    fprintf( '\r%s %i/%i', str, i, n )
end

