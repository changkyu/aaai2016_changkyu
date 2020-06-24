% ------------------------------------------------------------------------
% Generate data

clear; close all;

% Choose random seed: optional setting to reproduce numbers.
s = RandStream('mt19937ar','Seed',1);
RandStream.setGlobalStream(s);
reset(s,0);

% 20 dimensional, 500 samples searching, 5 dimensional subspace
VarW = 1;
VarMU = 1;
VarX = rand(1);

mW = 0;
mMU = 0;
mX = 0;

N = 1000;
D = 20;
M = 5;

% Z  (dim=M)   comes from N(0,I) 
Z  = mvnrnd(zeros(N,M), eye(M)); % each row is a sample
% W  (dim=D,M) comes from N(mW,vW)
W = orth(mvnrnd(ones(D,M)*mW, eye(M)*VarW)) * diag(M:-1:1);
% MU (dim=D)   comes from N(mMU,vMU)
MU = repmat(mvnrnd(ones(1,D)*mMU, eye(D)*VarMU), [N,1]);
% E  (dim=D)   comes from N(mX,vX)
%E = mvnrnd(ones(N,D)*mX, eye(D)*VarX);
E = ones(N,D)*mX + VarX*randn(N, D);
% model
X = W * Z' + MU' + E';

% ETA: Learning rate
ETA = 10;

% Node assignment to each sample
NV = 50;
V = get_sample_assign(NV, N);

% E: Adjacency graph topology (1:complete, 2:ring, 3:star, 4:chain)
GRAPHS = get_adj_graph(NV);
E = GRAPHS{2};

% For structure from motion applications only, better to enforce zero mean
ZeroMean = false;

% Convergence thresholds (typically you don't have to change these)
THRESHc = 1e-3;     % threshold (objective) for PPCA
THRESHd = 1e-3;     % threshold (objective) for D-PPCA
THRESHda = 1e-4;    % threshold (absolute)  for D-PPCA
THRESHdr = 1e-2;    % threshold (relative)  for D-PPCA

% Show progress per this number of iterations
objfreq_c = 1;      % show objective per this iterations (PPCA)
objfreq_d = 1;      % show objective per this iterations (D-PPCA)

% centralized model initialization: 
%   a random projection, zero mean and a small variance
m_init_c = get_init_value(X, M, 'ModelType', 'cppca');

% distributed model initialization:
%   initializes with a local PPCA result
m_init_d = get_init_value(X, M, 'ModelType', 'dppca', ...
    'NumberOfNodes', NV, 'SampleAssignVec', V);

% Run PPCA
disp('* PPCA (EM) *');
cm1 = cppca_em(X, M, ...
    'InitModel', m_init_c, 'Threshold', THRESHc, 'ShowObjPer', objfreq_c, ...
    'ZeroMean', ZeroMean, 'MaxIter', 1000);

% Run D-PPCA
disp('* D-PPCA *');
cm2 = dppca(X, M, V, E, 'Eta', ETA, ...
    'InitModel', m_init_d, 'Threshold', THRESHd, 'ShowObjPer', objfreq_d, ...
    'ThresholdA', THRESHda, 'ThresholdR', THRESHdr, ...
    'ZeroMean', ZeroMean, 'MaxIter', 1000, ...
    'UseResidual', true);

% Do D-PPCA with penalty parameter change (He, et al. 2000)
disp('* D-PPCA (varying penalty) *');
cm3 = dppca(X, M, V, E, 'Eta', ETA, ...
    'InitModel', m_init_d, 'Threshold', THRESHd, 'ShowObjPer', objfreq_d, ...
    'ThresholdA', THRESHda, 'ThresholdR', THRESHdr, ...
    'ZeroMean', ZeroMean, 'MaxIter', 1000, 'VaryETA', true, ...
    'VaryETA_mu', 0.1, 'VaryETA_tau', 1, 'VaryETA_max', 50, ...
    'UseResidual', true);

% Do D-PPCA with acceleration (Goldstein, et al. 2014)
disp('* D-PPCA (acceleration) *');
cm4 = dppca(X, M, V, E, 'Eta', ETA, ...
    'InitModel', m_init_d, 'Threshold', THRESHd, 'ShowObjPer', objfreq_d, ...
    'ThresholdA', THRESHda, 'ThresholdR', THRESHdr, ...
    'ZeroMean', ZeroMean, 'MaxIter', 1000, 'FwRes', true, ...
    'UseResidual', true);

% Plot
subplot(1,3,1);
if max([cm1.eITER, cm2.eITER, cm3.eITER]) < 100
    toPlot = cm1.objArray(1:cm1.eITER);
    toPlot = abs((toPlot - [realmax; toPlot(1:cm1.eITER-1)]) ./ [realmax; toPlot(1:cm1.eITER-1)]);
    semilogy(toPlot, 'gs-'); hold on;
    toPlot = cm2.objArray(1:cm2.eITER,NV+1);
    toPlot = abs((toPlot - [realmax; toPlot(1:cm2.eITER-1)]) ./ [realmax; toPlot(1:cm2.eITER-1)]);
    semilogy(toPlot, 'rx-'); 
    toPlot = cm3.objArray(1:cm3.eITER,NV+1);
    toPlot = abs((toPlot - [realmax; toPlot(1:cm3.eITER-1)]) ./ [realmax; toPlot(1:cm3.eITER-1)]);
    semilogy(toPlot, 'b-');
    toPlot = cm4.objArray(1:cm4.eITER,NV+1);
    toPlot = abs((toPlot - [realmax; toPlot(1:cm4.eITER-1)]) ./ [realmax; toPlot(1:cm4.eITER-1)]);
    semilogy(toPlot, 'mo-');
    %grid on;
else
    toPlot = cm1.objArray(1:cm1.eITER);
    toPlot = abs((toPlot - [realmax; toPlot(1:cm1.eITER-1)]) ./ [realmax; toPlot(1:cm1.eITER-1)]);
    loglog(toPlot, 'gs-'); hold on;
    toPlot = cm2.objArray(1:cm2.eITER,NV+1);
    toPlot = abs((toPlot - [realmax; toPlot(1:cm2.eITER-1)]) ./ [realmax; toPlot(1:cm2.eITER-1)]);
    loglog(toPlot, 'rx-'); 
    toPlot = cm3.objArray(1:cm3.eITER,NV+1);
    toPlot = abs((toPlot - [realmax; toPlot(1:cm3.eITER-1)]) ./ [realmax; toPlot(1:cm3.eITER-1)]);
    loglog(toPlot, 'b-');
    toPlot = cm4.objArray(1:cm4.eITER,NV+1);
    toPlot = abs((toPlot - [realmax; toPlot(1:cm4.eITER-1)]) ./ [realmax; toPlot(1:cm4.eITER-1)]);
    loglog(toPlot, 'mo-');
    %grid on;
end
legend({'Centralized', 'D-PPCA', 'D-PPCA (VP)', 'D-PPCA (F-ADMM)'});
ylabel('absolute relative change in objective');
xlabel('iterations');

subplot(1,3,2);
if max([cm1.eITER, cm2.eITER, cm3.eITER]) < 100
    plot(cm2.rArray(1:cm2.eITER), 'rx-'); hold on;  
    plot(cm2.rtArray(1:cm2.eITER), 'm--'); 
    plot(cm3.rArray(1:cm3.eITER), 'b-');
    plot(cm3.rtArray(1:cm3.eITER), 'c--');
    plot(cm4.rArray(1:cm4.eITER), 'g-');
    plot(cm4.rtArray(1:cm4.eITER), 'k--');
else
    semilogx(cm2.rArray(1:cm2.eITER), 'rx-'); hold on;  
    semilogx(cm2.rtArray(1:cm2.eITER), 'm--'); 
    semilogx(cm3.rArray(1:cm3.eITER), 'b-');
    semilogx(cm3.rtArray(1:cm3.eITER), 'c--');
    semilogx(cm4.rArray(1:cm4.eITER), 'g-');
    semilogx(cm4.rtArray(1:cm4.eITER), 'k--');
end
legend({'D-PPCA', 'D-PPCA (eps)', 'D-PPCA (VP)', 'D-PPCA (VP) (eps)', 'D-PPCA (F-ADMM)', 'D-PPCA (F-ADMM) (eps)'});
ylabel('primal residual');
xlabel('iterations');

subplot(1,3,3);
if max([cm1.eITER, cm2.eITER, cm3.eITER]) < 100
    plot(cm2.sArray(1:cm2.eITER), 'rx-'); hold on;  
    plot(cm2.stArray(1:cm2.eITER), 'm--'); 
    plot(cm3.sArray(1:cm3.eITER), 'b-');
    plot(cm3.stArray(1:cm3.eITER), 'c--');
    plot(cm4.sArray(1:cm4.eITER), 'go-');
    plot(cm4.stArray(1:cm4.eITER), 'k--');
else
    semilogx(cm2.sArray(1:cm2.eITER), 'rx-'); hold on;  
    semilogx(cm2.stArray(1:cm2.eITER), 'm--'); 
    semilogx(cm3.sArray(1:cm3.eITER), 'b-');
    semilogx(cm3.stArray(1:cm3.eITER), 'c--');
    semilogx(cm4.sArray(1:cm4.eITER), 'go-');
    semilogx(cm4.stArray(1:cm4.eITER), 'k--');
end
legend({'D-PPCA', 'D-PPCA (eps)', 'D-PPCA (VP)', 'D-PPCA (VP) (eps)', 'D-PPCA (F-ADMM)', 'D-PPCA (F-ADMM) (eps)'});
ylabel('dual residual');
xlabel('iterations');
