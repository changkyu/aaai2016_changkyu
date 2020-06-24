%--------------------------------------------------------------------------
% Variance correction for MFVB (Giordano and Broderick, ICML 2015)
%
% Implemented/Modified
%  by     Sejong Yoon (sjyoon@cs.rutgers.edu)
%  on     2014.04.02 (last modified on 2015/04/03)
%--------------------------------------------------------------------------
function [vW, vMU, vZ] = mfvb_correct(X, mW, mZ, mMU, vW, vZ, vMU, ...
    PX, PW, PMU)

[D, N] = size(X);
[~, M] = size(mW);

% Identities
I_W = speye(D*M);
I_Z = speye(M*N);
I_Mu = speye(D);
I_MuZ = speye(D + M*N);
I_WMu = speye(D*M + D);
I_WZ = speye(D*M + M*N);

% Variational covariance
V_W = sparse(diag(reshape(vW, [D*M, 1])));
V_Z = speye(M*N);
for n = 1 : N
    lb = M*(n-1)+1;
    ub = M*n;
    V_Z(lb:ub, lb:ub) = diag(diag(vZ(:,:,n)));
end
V_Mu = sparse(diag(vMU));
V_MuZ = sparse([ ...
    V_Mu                spalloc(D,M*N,0); ...
    spalloc(M*N,D,0),   V_Z ...
    ]);
V_WMu = sparse([ ...
    V_W                 spalloc(D*M,D,0); ...
    spalloc(D,D*M,0),   V_Mu ...
    ]);
V_WZ = sparse([ ...
    V_W                 spalloc(D*M,M*N,0); ...
    spalloc(M*N,D*M,0), V_Z ...
    ]);

% Hessians (symmetric)
% Eq. 5
H_W = ones(D*M, D*M);
BigZ = zeros(M,M); % Eq. 6
for n = 1 : N
    BigZ = BigZ + ( mZ(:,n) * mZ(:,n)' + inv(vZ(:,:,n)) );
end
for d1 = 1 : D
    for m1 = 1 : M
        for d2 = 1 : D
            for m2 = 1 : M
                % get index
                d1m1 = (m1 - 1) * D + d1;
                d2m2 = (m2 - 1) * D + d2;
                if d1 == d2 && m1 == m2
                    H_W(d1m1, d2m2) = -PX * BigZ(m1,m2) -PW(d1,m1);
                else
                    H_W(d1m1, d2m2) = -PX * BigZ(m1,m2);
                end
            end
        end
    end
end
H_W = sparse(H_W);
% Eq. 7
H_WMu = zeros(D*M, D);
for d1 = 1 : D
    for m = 1 : M
        for d2 = 1 : D
            if d1 == d2
                % get index
                d1m = (m - 1) * D + d1;
                H_WMu(d1m, d2) = -PX * sum(mZ(m,:));
            end
        end
    end
end
H_WMu = sparse(H_WMu);
% Eq. 8
H_WZ = zeros(D*M, M*N);
for d = 1 : D
    for m1 = 1 : M
        for m2 = 1 : M
            for n = 1 : N
                if m1 == m2
                    % get index
                    dm1 = (m1 - 1) * D + d;
                    m2n = (n - 1) * M + m2;
                    H_WZ(dm1, m2n) = -PX * sum(mW(d,:) .* mZ(:,n)') ...
                        + PX * (X(d,n) - mMU(d));
                end
            end
        end
    end
end
H_WZ = sparse(H_WZ);
% Eq. 9
H_Mu = speye(D);
H_Mu = H_Mu .* (- N * PX - PMU);
% Eq. 10 (just concatenating mW on the right)
H_MuZ = sparse(repmat(mW, [1, N]));
% Eq. 11
H_Z = zeros(M*N, M*N);
for m1 = 1 : M
    for n1 = 1 : N
        for m2 = 1 : M
            for n2 = 1 : N
                if n1 == n2
                    % get index
                    m1n1 = (n1 - 1) * M + m1;
                    m2n2 = (n2 - 1) * M + m2;
                    if m1 == m2
                        H_Z(m1n1,m2n2) = -PX * sum(mW(:,m1).*mW(:,m2)) - 1;
                    else
                        H_Z(m1n1,m2n2) = -PX * sum(mW(:,m1).*mW(:,m2));
                    end
                end
            end
        end
    end
end
H_Z = sparse(H_Z);
% These are transpose of the blocks computed above
H_MuW = H_WMu';
H_ZW = H_WZ';
H_ZMu = H_MuZ';

% Now we need to correct!
% Eq. 12
SIG_W = ( ...
    I_W - V_W * H_W ...
    - ( V_W * [H_WMu H_WZ] ) / ( I_MuZ - V_MuZ * [H_Mu H_MuZ; H_ZMu H_Z] ) ...
    * V_MuZ * [H_MuW; H_ZW] ...
    ) \ V_W;

% Eq. 13
SIG_Z = ( ...
    I_Z - V_Z * H_Z ...
    - ( V_Z * [H_ZW H_ZMu] ) / ( I_WMu - V_WMu * [H_W H_WMu; H_MuW H_Mu] ) ...
    * V_WMu * [H_WZ; H_MuZ] ...
    ) \ V_Z;

% Eq. 14
SIG_Mu = ( ...
    I_Mu - V_Mu * H_Mu ...
    - ( V_Mu * [H_MuW H_MuZ] ) / ( I_WZ - V_WZ * [H_W H_WZ; H_ZW H_Z] ) ...
    * V_WZ * [H_WMu; H_ZMu] ...
    ) \ V_Mu;

% Now update the variational covariances
vW = full(reshape(diag(SIG_W), [D, M]));
vZ = zeros(M,M,N);
for n = 1 : N
    lb = M*(n-1)+1;
    ub = M*n;
    vZ(:,:,n) = full(diag(diag(SIG_Z(lb:ub, lb:ub))));
end
vMU = full(diag(SIG_Mu));

end