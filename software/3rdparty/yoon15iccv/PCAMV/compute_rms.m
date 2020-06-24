%  COMPUTE_RMS - RMS error and reconstruction error matrix
%
%  [rms,errMx] = COMPUTE_RMS( X, A, S, M, ndata ) computes a matrix
%  errMx of reconstruction errors (X - A*S) and the rms reconstruction
%  error.
%
%  [rms,errMx] = COMPUTE_RMS( ..., numCPU ) also specifies the number
%  numCPU of CPUs used for parallel computing (default 1).
%
%  See also ERRPCA_PT.CPP

%  This software is provided "as is", without warranty of any kind.
%  Alexander Ilin, Tapani Raiko

function [rms,errMx] = compute_rms( X, A, S, M, ndata, numCPU )

if isempty(X)
    errMx = []; rms = NaN; return
end

if issparse(X)
    if nargin < 6, numCPU = 1; end
    
    %Linux:
    errMx = errpca_pt( X, A, S, numCPU );
else
    errMx = (X - A*S).*M;
end

rms = full(sqrt(sum(sum(errMx.^2))/ndata));

