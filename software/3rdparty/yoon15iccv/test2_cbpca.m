%--------------------------------------------------------------------------
% Test on our data

clear; close all;

% Path
addpath(genpath('PCAMV'));

% Choose random seed: optional setting to reproduce numbers.
s = RandStream('mt19937ar','Seed',1);
RandStream.setGlobalStream(s);
reset(s,0);

% 50 dimensional, 300 samples searching, 5 dimensional subspace
VarW = rand(1);
VarMU = rand(1);
VarX = rand(1);

mW = 10;
mMU = 5;
mX = 0;

N = 1000;
D = 50;
M = 10;

% Z  (dim=M)   comes from N(0,I) 
Z = mvnrnd(zeros(N,M), eye(M)); % each row is a sample
% W  (dim=D,M) comes from N(mW,vW)
W = mvnrnd(ones(D,M)*mW, eye(M)*VarW);
% MU (dim=D)   comes from N(mMU,vMU)
MU = repmat(mvnrnd(ones(1,D)*mMU, eye(D)*VarMU), [N,1]);
% E  (dim=D)   comes from N(mX,vX)
%E = mvnrnd(ones(N,D)*mX, eye(D)*VarX);
E = ones(N,D)*mX + VarX*randn(N, D);
% model
X = W * Z' + MU' + E';

% PPCA (Ilin)
disp('PPCA (Ilin)');
opts = struct( ...
    'algorithm', 'ppca',... % algorithm
    'maxiters', 1000,... % maximum iterations
    'cfstop', [ 1 0 1e-5 ], ... % stop if cost not change over some iters.
    ... % [#iterations_to_monitor, absolute_error, relative_error]
    'minangle', 0, ... % keep going regardless of subspace angle change
    'uniquesv', 0 ... % compute only unique covariance matrices of S(:,j)
);
output = pca_full( X, M, opts );
cm1 = output2model(output);

% PPCA (ours)
disp('PPCA (ours)');
cm2 = cppca_em( X, M );

% VBPCA (Ilin)
disp('VBPCA (Ilin)');
opts.algorithm = 'vb';
% The prior is not updated at the begining of learning to avoid killing sources 
opts.niter_broadprior = 5; 
output = pca_full( X, M, opts );
cm3 = output2model(output);

% VBPCA (Ours)
disp('VBPCA (Ours)');
cm4 = cbpca( X, M );

% Z estimate
disp('* Z estimates (root mean square error):');
rms_res(1) = rms(rms(Z - cm1.EZ'));
rms_res(2) = rms(rms(Z - cm2.EZ'));
rms_res(3) = rms(rms(Z - cm3.mZ'));
rms_res(4) = rms(rms(Z - cm4.mZ'));
fprintf('            1e-123456789012345\n');
fprintf('PPCA-Ilin  : %f\n', rms_res(1));
fprintf('PPCA-ours  : %f\n', rms_res(2));
fprintf('VBPCA-Ilin : %f\n', rms_res(3));
fprintf('VBPCA-ours : %f\n', rms_res(4));

% X reconstruction
disp('* X reconstruction (root mean square error):');
fprintf('            1e-123456789012345\n');
fprintf('PPCA-Ilin  : %.15f\n', rms(rms(X - (cm1.W * cm1.EZ + repmat(cm1.MU, [1, N])))));
fprintf('PPCA-ours  : %.15f\n', rms(rms(X - (cm2.W * cm2.EZ + repmat(cm2.MU, [1, N])))));
fprintf('VBPCA-Ilin : %.15f\n', rms(rms(X - (cm3.mW * cm3.mZ + repmat(cm3.mMU, [1, N])))));
fprintf('VBPCA-ours : %.15f\n', rms(rms(X - (cm4.mW * cm4.mZ + repmat(cm4.mMU, [1, N])))));

% W estimate
disp('* W estimates (root mean square error):');
fprintf('            1e-123456789012345\n');
fprintf('PPCA-Ilin  : %.15f\n', rms(rms(W - cm1.W)));
fprintf('PPCA-ours  : %.15f\n', rms(rms(W - cm2.W)));
fprintf('VBPCA-Ilin : %.15f\n', rms(rms(W - cm3.mW)));
fprintf('VBPCA-ours : %.15f\n', rms(rms(W - cm4.mW)));

% W estimate
disp('* W estimates (subspace angle):');
fprintf('            1e-123456789012345\n');
fprintf('PPCA-Ilin  : %.15f\n', subspace(W, cm1.W));
fprintf('PPCA-ours  : %.15f\n', subspace(W, cm2.W));
fprintf('VBPCA-Ilin : %.15f\n', subspace(W, cm3.mW));
fprintf('VBPCA-ours : %.15f\n', subspace(W, cm4.mW));

% MU estimate
disp('* MU estimates (root mean square error):');
fprintf('            1e-123456789012345\n');
fprintf('PPCA-Ilin  : %.15f\n', rms(MU(1,:)' - cm1.MU));
fprintf('PPCA-ours  : %.15f\n', rms(MU(1,:)' - cm2.MU));
fprintf('VBPCA-Ilin : %.15f\n', rms(MU(1,:)' - cm3.mMU));
fprintf('VBPCA-ours : %.15f\n', rms(MU(1,:)' - cm4.mMU));

% Variance estimate
disp('* Variance estimates (absolute error):');
fprintf('GT variance: %f\n', VarX);
fprintf('            1e-123456789012345\n');
fprintf('PPCA-Ilin  : %.15f\n', abs(VarX - cm1.VAR));
fprintf('PPCA-ours  : %.15f\n', abs(VarX - cm2.VAR));
fprintf('VBPCA-Ilin : %.15f\n', abs(VarX - cm3.PX));
fprintf('VBPCA-ours : %.15f\n', abs(VarX - cm4.PX));

% PW and PMU are point estimates on the hyperparameters (VarW, VarMU)
% You can check them by cm3.PW, cm3.PMU for their implementation and 
% cm4.PW and cm4.PMU for our implementation.