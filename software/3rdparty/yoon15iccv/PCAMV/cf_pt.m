%  CF_PT - Cost for PCA model with maximum a posteriori estimation
%
%  Compute the cost function for the PCA model with maximum a
%  posteriori estimation of the parameters. The function also returns
%  the matrix of reconstruction errors errMx and the rms error.
%
%  See also PCA_PT

%  This software is provided "as is", without warranty of any kind.
%  Alexander Ilin, Tapani Raiko

function [cost,errMx,rms] = ...
    cf_pt( X, A, S, Mu, V, Va, Vmu, M, ndata, numCPU )

if issparse(X)
    % X is a sparse matrix with only observed values
    errMx = errpca_pt( X, A, S, numCPU );
else
    % Missing values are marked as NaNs
    errMx = (X - A*S).*M;
end

[n1,n2] = size(X);
ncomp = size(A,2);

err2 = full(sum(sum(errMx.^2)));
rms = sqrt(err2/ndata);

cost = 0.5/V*err2 + ndata*0.5*log(2*pi*V);

if all(~isinf(Va))
    % Prior for A, Mu and S
    if Vmu ~= 0
        cost = cost + 0.5/Vmu*sum(Mu.^2) + n1*0.5*log(2*pi*Vmu);
    end

    cost = cost + 0.5*sum(sum(A.^2,1)./Va) + n1*0.5*sum(log(2*pi*Va));
    
    cost = cost + 0.5*sum(sum(S.^2)) + n2*ncomp*0.5*log(2*pi);
    
end

