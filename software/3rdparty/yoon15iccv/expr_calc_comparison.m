function [res] = expr_calc_comparison(m_cppca, m_cbpca, m_dppca, m_dbpca, ...
    X, W, Z, MU, mX, mW, mMU, VarX, VarW, VarMU)

cm1 = m_cppca;
cm2 = m_cbpca;
cm3 = m_dppca;
cm4 = m_dbpca;

[D, N] = size(X);   % data dimension and number of samples
M = size(W, 2);     % subspace dimension
NV = length(cm3.W); % number of nodes for distributed models

% OUTPUT
fprintf('##### RESULT (Root Mean Square vs. Ground Truth) #####\n');
fprintf('Examples : %d / Dimension: %d / Subspace Dimension: %d\n', N, D, M);

%----------------------------------------------------------------------
fprintf('*** Ground Truth ***\n');
fprintf('bias     of X : %f\n', mX);
fprintf('bias     of W : %f\n', mW);
fprintf('bias     of MU: %f\n', mMU);
fprintf('variance of X : %f\n', VarX);
fprintf('variance of W : %f\n', VarW);
fprintf('variance of MU: %f\n', VarMU);

%----------------------------------------------------------------------
% Z estimate
disp('* Z estimates (root mean square error):');
rms_Z(1) = rms(rms(Z - cm1.EZ'));
rms_Z(2) = rms(rms(Z - cm2.mZ'));

rms_Z(3) = rms(rms(Z - cell2mat(cm3.EZ')'));
rms_Z(4) = rms(rms(Z - cell2mat(cm4.mZ')'));

fprintf('            1e-123456789012345\n');
fprintf('PPCA-ours  : %.15f\n', rms_Z(1));
fprintf('VBPCA-ours : %.15f\n', rms_Z(2));
fprintf('D-PPCA     : %.15f\n', rms_Z(3));
fprintf('D-VBPCA    : %.15f\n', rms_Z(4));

%----------------------------------------------------------------------
% X reconstruction
disp('* X reconstruction (root mean square error):');
rms_X(1) = calc_ppca_rms(X, cm1.W, cm1.EZ, cm1.MU);
rms_X(2) = calc_ppca_rms(X, cm2.mW, cm2.mZ, cm2.mMU);

rms_X(3) = calc_ppca_rms(X, cm3.W, cm3.EZ, cm3.MU);
rms_X(4) = calc_ppca_rms(X, cm4.mW, cm4.mZ, cm4.mMU);

fprintf('            1e-123456789012345\n');
fprintf('PPCA-ours  : %.15f\n', rms_X(1));
fprintf('VBPCA-ours : %.15f\n', rms_X(2));
fprintf('D-PPCA     : %.15f\n', rms_X(3));
fprintf('D-VBPCA    : %.15f\n', rms_X(4));

%----------------------------------------------------------------------
% W estimate
disp('* W estimates (root mean square error):');
rms_W(1) = rms(rms(W - cm1.W));
rms_W(2) = rms(rms(W - cm2.mW));

rms_W(3) = -Inf;
rms_W(4) = -Inf;
for i = 1 : NV
    t = rms(rms(W - cm3.W{i}));
    if rms_W(3) < t
        rms_W(3) = t;
    end
    t = rms(rms(W - cm4.mW{i}));
    if rms_W(4) < t
        rms_W(4) = t;
    end
end

fprintf('            1e-123456789012345\n');
fprintf('PPCA-ours  : %.15f\n', rms_W(1));
fprintf('VBPCA-ours : %.15f\n', rms_W(2));
fprintf('D-PPCA     : %.15f\n', rms_W(3));
fprintf('D-VBPCA    : %.15f\n', rms_W(4));

%----------------------------------------------------------------------
% W estimate
disp('* W estimates (subspace angle):');
rms_WA(1) = subspace(W, cm1.W);
rms_WA(2) = subspace(W, cm2.mW);

rms_WA(3) = -Inf;
rms_WA(4) = -Inf;
for i = 1 : NV
    t = subspace(W, cm3.W{i});
    if rms_WA(3) < t
        rms_WA(3) = t;
    end
    t = subspace(W, cm4.mW{i});
    if rms_WA(4) < t
        rms_WA(4) = t;
    end
end

fprintf('            1e-123456789012345\n');
fprintf('PPCA-ours  : %.15f\n', rms_WA(1));
fprintf('VBPCA-ours : %.15f\n', rms_WA(2));
fprintf('D-PPCA     : %.15f\n', rms_WA(3));
fprintf('D-VBPCA    : %.15f\n', rms_WA(4));

%----------------------------------------------------------------------
% MU estimate
disp('* MU estimates (root mean square error):');
rms_MU(1) = rms(MU(1,:)' - cm1.MU);
rms_MU(2) = rms(MU(1,:)' - cm2.mMU);

rms_MU(3) = -Inf;
rms_MU(4) = -Inf;
for i = 1 : NV
    t = rms(MU(1,:)' - cm3.MU{i});
    if rms_MU(3) < t
        rms_MU(3) = t;
    end
    t = rms(MU(1,:)' - cm4.mMU{i});
    if rms_MU(4) < t
        rms_MU(4) = t;
    end
end

fprintf('            1e-123456789012345\n');
fprintf('PPCA-ours  : %.15f\n', rms_MU(1));
fprintf('VBPCA-ours : %.15f\n', rms_MU(2));
fprintf('D-PPCA     : %.15f\n', rms_MU(3));
fprintf('D-VBPCA    : %.15f\n', rms_MU(4));

%----------------------------------------------------------------------
% Variance estimate
disp('* Variance estimates (absolute error):');
rms_VarX(1) = abs(VarX - cm1.VAR);
rms_VarX(2) = abs(VarX - cm2.PX);

rms_VarX(3) = -Inf;
rms_VarX(4) = -Inf;
for i = 1 : NV
    t = abs(VarX - cm3.VAR{i});
    if rms_VarX(3) < t
        rms_VarX(3) = t;
    end
    t = abs(VarX - cm4.PX{i});
    if rms_VarX(4) < t
        rms_VarX(4) = t;
    end
end

fprintf('GT variance: %f\n', VarX);
fprintf('            1e-123456789012345\n');
fprintf('PPCA-ours  : %.15f\n', rms_VarX(1));
fprintf('VBPCA-ours : %.15f\n', rms_VarX(2));
fprintf('D-PPCA     : %.15f\n', rms_VarX(3));
fprintf('D-VBPCA    : %.15f\n', rms_VarX(4));

res = structure(rms_Z, rms_X, rms_W, rms_WA, rms_MU, rms_VarX);

end
