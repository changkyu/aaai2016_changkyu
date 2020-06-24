%  CF_FULL - Cost for PCA model with unrestricted Gaussian posterior
%
%  cost = CF_FULL( X, A, S, Mu, V, cv, hp ) computes the cost function
%  for the PCA model with unrestricted Gaussian pdfs (full covariance
%  matrices) for A(i,:) and S(:,j). The input parameters are returned
%  by PCA_FULL.
%
%  See also PCA_FULL

%  The function is also called by PCA_FULL providing the parameters
%  that have already been computed.

%  This software is provided "as is", without warranty of any kind.
%  Alexander Ilin, Tapani Raiko

function [ cost, cost_x, cost_a, cost_mu, cost_s ] = ...
    cf_full( X, A, S, Mu, V, Av, Sv, Isv, Muv, Va, Vmu, M, sXv, ndata )

[n1,n2] = size(X);
ncomp = size(A,2);

if nargin == 7
    cv = Av;
    hp = Sv;
    
    Av = cv.A;
    Sv = cv.S;
    Isv = cv.Isv;
    Muv = cv.Mu;
    Va = hp.Va;
    Vmu = hp.Vmu;
    if issparse(X),  M = spones(X);
    else,            M = ~isnan(X); X(isnan(X)) = 0;
    end
    sXv = [];

    if ~isempty(Mu)
        X = X - repmat(Mu,1,n2).*M;
    end
end

if isempty(sXv)
    [IX,JX,data] = find(M);
    ndata = length(IX);
    [rms,errMx] = compute_rms( X, A, S, M, ndata );

    % TODO: This is inefficient?
    sXv = (rms^2)*ndata;
    if isempty(Isv)
        for r = 1:ndata
            sXv = sXv + A(IX(r),:) * Sv{JX(r)} * A(IX(r),:)';
            if ~isempty(Av)
                sXv = sXv + S(:,JX(r))' * Av{IX(r)} * S(:,JX(r)) ...
                      + sum( sum( Sv{JX(r)} .* Av{IX(r)} ) );
            end
        end
    else
        for r = 1:ndata
            sXv = sXv + A(IX(r),:) * Sv{Isv(JX(r))} * A(IX(r),:)';
            if ~isempty(Av)
                sXv = sXv + S(:,JX(r))' * Av{IX(r)} * S(:,JX(r)) ...
                      + sum( sum( Sv{Isv(JX(r))} .* Av{IX(r)} ) );
            end
        end
    end
    if ~isempty(Muv)
        sXv = sXv + sum(Muv(IX));
    end
end

if 0
    cost = 0;
    for j = 1:n2
        A_j = repmat(full(M(:,j)),1,ncomp) .* A;
        if isempty(Isv)
            S1 = S(:,j)*S(:,j)' + Sv{j};
            cost = cost - 0.5 * log(det(Sv{j})) ...
                   - ncomp/2*( 1 + log(2*pi) );
        else
            S1 = S(:,j)*S(:,j)' + Sv{Isv(j)};
            cost = cost - 0.5 * log(det(Sv{Isv(j)})) ...
                   - ncomp/2*( 1 + log(2*pi) );
        end
        
        cost = cost + ...
               sum(M(:,j),1)/2 * log(V) ...
               + 0.5 * trace( S1 ) ...
               + 0.5 / V * X(:,j)' * X(:,j) ...
               - 1 / V * S(:,j)' * A_j' * X(:,j) ...
               + 0.5 / V * trace( A_j'*A_j*S1 );
    end

    %This might be faster
    %cost2 = 0;
    %for j = 1:n2
    %    A_j = repmat(full(M(:,j)),1,ncomp) .* A;
    %    cost2 = cost2 + 0.5 / V * trace( A_j'*A_j*Sv{j} );
    %end

else
    use_prior = all(~isinf(Va));
    
    cost_x = 0.5 / V * sXv + 0.5*ndata*log( 2*pi*V );
    
    cost_mu = 0;
    cost_a = 0;

    if use_prior
        if ~isempty(Muv)
            cost_mu = 0.5/Vmu*sum(Mu.^2+Muv) ...
                      - 0.5*sum(log(Muv)) + n1/2*log(Vmu) - n1/2;
        elseif Vmu ~= 0
            cost_mu = 0.5/Vmu*sum(Mu.^2) + n1/2*log(2*pi*Vmu);
        end
        
        if ~isempty(Av)
            cost_a = 0.5*sum(sum(A.^2,1)./Va) ...
                     + n1/2*sum(log(Va)) - n1*ncomp/2;
            for i = 1:n1
                cost_a = cost_a ...
                         + 0.5*sum( diag(Av{i})' ./Va ) ...
                         - 0.5*log(det(Av{i}));
            end
        else
            cost_a = 0.5*sum(sum(A.^2,1)./Va) + n1/2*sum(log(2*pi*Va));
        end

    
    else % No prior for A and Mu

        %    cost = cost - 0.5*sum(sum(log(2*pi*Av))) - n1*ncomp/2;

        if ~isempty(Muv)
            cost_mu = -0.5*sum(log(2*pi*Muv)) - n1/2;
        end
        if ~isempty(Av)
            cost_a = -n1*ncomp/2*(1+log(2*pi));
            for i = 1:n1
                cost_a = cost_a - 0.5*log(det(Av{i}));
            end
        end
    
    end

    cost_s = 0.5*sum(sum(S.^2));
    if isempty(Isv)
        for j = 1:n2
            cost_s = cost_s + 0.5*trace(Sv{j}) - 0.5*log(det(Sv{j}));
        end
    else
        for j = 1:n2
            cost_s = cost_s + 0.5*trace(Sv{Isv(j)}) - ...
                     0.5*log(det(Sv{Isv(j)}));
        end
    end
    cost_s = cost_s - ncomp*n2/2;

    cost = cost_mu + cost_a + cost_x + cost_s;

end

if 0
    cost_x
    cost_a
    cost_mu
    cost_s
end
