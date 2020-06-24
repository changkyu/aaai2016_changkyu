%--------------------------------------------------------------------------
% Test on toy spiral subspace data (from Ilin and Raiko's code)

clear; close all;

% Path
addpath(genpath('PCAMV'));

% Choose random seed: optional setting to reproduce numbers.
s = RandStream('mt19937ar','Seed',1);
RandStream.setGlobalStream(s);
reset(s,0);

% Generate the matrix of inputs x and targets t.
N = 500;
D = 5;
M = 2;
VarX = 0.2;

W = orth(randn(D,M))*diag(M:-1:1);
T = 1:N;
Z = [ exp(-T/150).*cos( 2*pi*T/50 );
      exp(-T/150).*sin( 2*pi*T/50 ) ];
   
% Normalizing to zero mean and unit variance
Z = ( Z - repmat( mean(Z,2), 1, N ) );
Z = Z ./ repmat( sqrt( mean( Z.^2, 2 ) ), 1, N );
X = W * Z;
X = X + VarX * randn(D,N);

% Missing rate of data
MissRate = 1;       % missing rate in percentage
if MissRate > 0
    seq = randperm(D * N);
    seq = seq(1:floor(D * N * MissRate / 100));
    X(seq) = NaN;
end

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

% PPCA (Ours)
disp('PPCA (Ours)');
cm2 = cppca_em( X, M );

% VBPCA (Ilin)
disp('VBPCA (Ilin)');
opts.algorithm = 'vb';
output = pca_full( X, M, opts );
cm3 = output2model(output);

% VBPCA (Ours)
disp('VBPCA (Ours)');
cm4 = cbpca( X, M );

% Plot results
figure;
subplot(2,3,1); plot(Z(1,:),Z(2,:),'+'); 
axis equal; title('Z'); xlim([-3.5 3.5]); ylim([-3.5 3.5]);
subplot(2,3,2); plot(cm1.EZ(1,:),cm1.EZ(2,:),'o'); 
axis equal; title('PPCA (IR)'); xlim([-3.5 3.5]); ylim([-3.5 3.5]);
subplot(2,3,3); plot(cm2.EZ(1,:),cm2.EZ(2,:),'o'); 
axis equal; title('PPCA (ours)'); xlim([-3.5 3.5]); ylim([-3.5 3.5]);
subplot(2,3,5); plot(cm3.mZ(1,:),cm3.mZ(2,:),'o'); 
axis equal; title('VBPCA (IR)'); xlim([-3.5 3.5]); ylim([-3.5 3.5]);
subplot(2,3,6); plot(cm4.mZ(1,:),cm4.mZ(2,:),'o'); 
axis equal; title('VBPCA (ours)'); xlim([-3.5 3.5]); ylim([-3.5 3.5]);

figure;
for d = 1 : D
    subplot(D,1,d); 
    
    Xrpj1 = bsxfun(@plus, cm1.W * cm1.EZ, cm1.MU);
    Xrpj2 = bsxfun(@plus, cm2.W * cm2.EZ, cm2.MU);
    Xrpj3 = bsxfun(@plus, cm3.mW * cm3.mZ, cm3.mMU);
    Xrpj4 = bsxfun(@plus, cm4.mW * cm4.mZ, cm4.mMU);
    data = [X(d,:); Xrpj1(d,:); Xrpj2(d,:); Xrpj3(d,:); Xrpj4(d,:)];
    
    plot(data');
    ylim([min(min(data)) max(max(data))]);
    legend('X', 'PPCA (IR)', 'PPCA (ours)', 'VBPCA (IR)', 'VBPCA (ours)', ...
        'location', 'eastoutside');
end

% Subspace angle of W
disp('* Subspace angle of W (vs. GT):');
fprintf('            1e-123456789012345\n');
fprintf('PPCA-Ilin  : %.15f\n', subspace(W,cm1.W));
fprintf('PPCA-ours  : %.15f\n', subspace(W,cm2.W));
fprintf('VBPCA-Ilin : %.15f\n', subspace(W,cm3.mW));
fprintf('VBPCA-ours : %.15f\n', subspace(W,cm4.mW));

% Mean estimate
disp('* Mean estimates (MU) :');
disp('GT        PPCA-Ilin PPCA-ours BPCA-Ilin BPCA-ours');
disp([mean(X,2) cm1.MU cm2.MU cm3.mMU cm4.mMU]);

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