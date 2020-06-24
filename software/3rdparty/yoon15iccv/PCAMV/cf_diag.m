%  CF_DIAG - Cost for PCA model with fully factorial posterior
%
%  cost = CF_DIAG( X, A, S, Mu, V, cv, hp ) computes the cost function
%  for the PCA model with fully factorial Gaussian pdf (diagonal
%  covariance matrix) for approximating the posterior distributions of
%  A and S. The input parameters are returned by PCA_DIAG.
%
%  [cost,errMx,rms] = CF_DIAG( ) also returns the matrix of
%  reconstruction errors errMx and the rms error.
%
%  See also PCA_DIAG

%  The function is also called by PCA_DIAG providing the parameters
%  that have already been computed.

%  This software is provided "as is", without warranty of any kind.
%  Alexander Ilin, Tapani Raiko

function [cost,errMx,rms,vN] = ...
    cf_diag( X, A, S, Mu, V, Av, Sv, Muv, Va, Vmu,...
             M, Nobs_i, ndata, numCPU )

[n1,n2] = size(X);
ncomp = size(A,2);

if nargin == 7
    cv = Av;
    hp = Sv;
    
    Av = cv.A;
    Sv = cv.S;
    Muv = cv.Mu;
    Va = hp.Va;
    Vmu = hp.Vmu;
    if issparse(X),  M = spones(X);
    else,            M = ~isnan(X); X(isnan(X)) = 0;
    end
    Nobs_i = full(sum(M,2));
    ndata = sum(Nobs_i);
    numCPU = 1;

    if ~isempty(Mu)
        X = X - repmat(Mu,1,n2).*M;
    end

end

if issparse(X)
    % X is a sparse matrix with only observed values
    [errMx,varcost] = errpca_diag( X, A, S, Sv, Av, numCPU );
else
    % Missing values are marked as NaNs
    errMx = (X - A*S).*M;
    if isempty(Av)
        varcost = ( (A.^2)*Sv ).*M;
    else
        varcost = ( Av*S.^2 + (A.^2)*Sv + Av*Sv ).*M;
    end
    varcost = full(sum(sum(varcost)));
end

err2 = sum(full(sum(errMx.^2,2)));
vN = err2 + varcost;
if ~isempty(Muv)
    vN = vN + sum(Nobs_i.*Muv);
end

rms = sqrt(err2/ndata);

cost_x = 0.5/V*vN + ndata*0.5*log(2*pi*V);

if all(~isinf(Va))
    % Prior for A and Mu
    if ~isempty(Muv)
        % VB
        cost_mu = 0.5/Vmu*sum(Mu.^2+Muv) ...
                  - 0.5*sum(log(Muv)) + n1/2*log(Vmu) - n1/2;
    elseif Vmu ~= 0
        cost_mu = 0.5/Vmu*sum(Mu.^2) + n1/2*log(2*pi*Vmu);
    end
    
    if ~isempty(Av)
        cost_a = 0.5*sum(sum(A.^2+Av,1)./Va) ...
                 - 0.5*sum(sum(log(Av))) + n1/2*sum(log(Va)) - n1*ncomp/2;
    else
        cost_a = 0.5*sum(sum(A.^2,1)./Va) + n1/2*sum(log(2*pi*Va));
        %cost = cost + 0.5*sum(sum(A.^2,1)./Va) + ...
        %       n1/2*sum(log(Va)) + n1*ncomp/2*log(2*pi);
    end
    
else
    % No prior for A and Mu
    if ~isempty(Muv)
        %cost = cost - 0.5*sum(log(2*pi*Muv)) - n1/2;
        cost_mu = - 0.5*sum(log(Muv)) - n1/2*(1+log(2*pi));
    else
        cost_mu = 0;
    end
    if ~isempty(Av)
        %cost = cost - 0.5*sum(sum(log(2*pi*Av))) - n1*ncomp/2;
        cost_a = - 0.5*sum(sum(log(Av))) - n1*ncomp/2*(1+log(2*pi));
    else
        cost_a = 0;
    end
end

cost_s = 0.5*sum(sum(S.^2+Sv)) - 0.5*sum(sum(log(Sv))) - n2*ncomp/2;

cost = cost_x + cost_mu + cost_a + cost_s;

if 0
    cost_x
    cost_a
    cost_mu
    cost_s
end
