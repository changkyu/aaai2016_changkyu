%  PCAIMPUT - Imputation algorithm for PCA with missing values
%
%  [ A, S, Mu, LC ] = PCAIMPUT( X, N ) finds the PCA decomposition X =
%  AS + [Mu .. Mu] for the given data matrix X and number of principal
%  components N using repetitive imputation. Matrix X can be either
%  sparse with only observed values included (note that observed zeros
%  should be replaced with eps) or a full matrix with missing values
%  replaced by NaNs.
%
%  LC is a structure with learning curves (rms training and probing
%  errors, time).
%
%  Optional parameter/value pairs with {default values}:
%   init       - Structure with initial parameter values {[]}
%   maxiters   - Maximum number of iterations {1000}
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
%   filename   - Name of the file for auto-save {'pcaimp_autosave'}
%
%  OUT = PCAIMPUT( X, N ) returns all the outputs in a single
%  structure variable OUT.

%  This software is provided "as is", without warranty of any kind.
%  Alexander Ilin, Tapani Raiko

function [ A, S, Mu, lc ] = pcaimput( X, ncomp, varargin )

opts = struct( ...
    'init',          [],...
    'maxiters',      1000,...
    'minangle',      0,...
    'bias',          1,...
    'autosave',      3600,...
    'filename',      'pcaimp_autosave',...
    'earlystop',     0,...
    'rmsstop',       [],... % [] means no rms stop criteria
    'verbose',       1,...
    'xprobe',        [],...
    'display',       0 );

[ opts, errmsg, wrnmsg ] = argschk( opts, varargin{:} );
if ~isempty(errmsg), error( errmsg ), end
if ~isempty(wrnmsg), warning( wrnmsg ), end
Xprobe = opts.xprobe;
switch opts.verbose
case {0,1}
    eigsopts.disp = 0;
case 2
    eigsopts.disp = 1;
end

[n1x,n2x] = size(X);
[ X, Xprobe, Ir, Ic, opts.init ] = rmempty( X, Xprobe,...
                                            opts.init, opts.verbose );
[n1,n2] = size(X);

if issparse(X)
    % X is a sparse matrix with only observed values
    M = spones(X);    Mprobe = spones(Xprobe);
    X = full(X);      Xprobe = full(Xprobe);
else
    % Missing values are marked as NaNs
    M = ~isnan(X);    Mprobe = ~isnan(Xprobe);

    X(X==0) = eps;    Xprobe(Xprobe==0) = eps;
    X(isnan(X)) = 0;  Xprobe(isnan(Xprobe)) = 0;
end
Mmis = ~M;
Nobs_i = sum(M,2);
ndata = sum(Nobs_i);

nprobe = nnz(Mprobe);
if nprobe == 0
    Xprobe = [];
    opts.earlystop = 0;
    prms = NaN;
end

if ndata == n1*n2
    opts.maxiters = 1;
end

Xorig = X;
lc.rms = NaN; lc.prms = NaN; lc.time = 0;

[ A, S, Mu ] = InitParms( opts.init );

if isempty(Mu)
    if ~isempty(A) && ~isempty(S)
        X = zeros(size(X));
    else
        % Replace missing values with the row-wise mean
        Mu = sum(X,2)./Nobs_i;
        X = repmat(Mu,1,n2);
    end
end
if ~isempty(A) && ~isempty(S)
    X = X + A*S;
end
X = Xorig + X.*Mmis;

time_start = clock;
time_autosave = time_start;
tic
for iter = 1:opts.maxiters

    % Remove the mean
    if opts.bias
        Mu = mean(X,2);
        X = X - repmat(Mu,1,n2);
        Xprobe = X - repmat(Mu,1,n2);
    end
    
    % Compute PCA
    C = X*X'/n2;
    [A,D] = eigs( C, ncomp, 'LM', eigsopts );
    S = A'*X;

    X = repmat(Mu,1,n2) + A*S;
    rms = sqrt(sum(sum((Xorig-X.*M).^2))/ndata);
    if nprobe > 0
        prms = sqrt(sum(sum((Xprobe-X.*Mprobe).^2))/nprobe);
    end

    % Replace missing values with mu + A*S
    if iter < opts.maxiters
        X = Xorig + X.*Mmis;
    end
    
    t = toc;
    lc.rms = [ lc.rms rms ]; lc.prms = [ lc.prms prms ];
    lc.time = [ lc.time t ];

    if iter > 1
        angleA = subspace(A,Aold);
    else
        angleA = NaN;
    end
    PrintStep( opts.verbose, lc, angleA )

    convmsg = converg_check( opts, lc, angleA );
    if ~isempty(convmsg)
        if opts.verbose, fprintf( '%s', convmsg ), end
        break
    end
    Aold = A;
    
    time = clock;
    if etime(time,time_autosave) > opts.autosave
        time_autosave = time;
        save( opts.filename, 'A', 'S', 'Mu', 'lc',...
              'Ir', 'Ic', 'n1x', 'n2x', 'n1', 'n2' )
    end
end

D = diag(D);
S = diag(1./sqrt(D))*S;
A = A*diag(sqrt(D));

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
function [ A, S, Mu ] = InitParms( init )

A = []; S = []; Mu = [];

if isstruct(init)
    if isfield( init, 'A' )
        A = init.A;
    end
    if isfield( init, 'Mu' )
        Mu = init.Mu;
    end
    if isfield( init, 'S' )
        S = init.S;
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

