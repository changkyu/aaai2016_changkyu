function model = output2model(output)
% Convert PPCA / VBPCA output of Ilin and Raiko (2010) into our format

% output structure example for D = 5, M = 2
%       A: [5x2 double]
%       S: [2x500 double]
%      Mu: [5x1 double]
%       V: 0.0385
%      Va: [1000 1000]
%     Vmu: 1000
%      Av: {[2x2 double]  [2x2 double]  [2x2 double]  [2x2 double]  [2x2 double]}
%      Sv: {1x500 cell}
%     Isv: []
%     Muv: [5x1 double]
%      lc: [1x1 struct]
[D, M] = size(output.A);
[~, N] = size(output.S);

% PPCA? VBPCA?
if isempty(output.Av)
%  W        : D x M projection matrix
%  MU       : D x 1 vector sample means
%  VAR      : Scalar estimated variance
%  EZ       : M x N matrix mean of N latent vectors
%  EZZt     : M x M x N cube for covariance of N latent vectors
    model.W = output.A;
    model.MU = output.Mu;
    model.VAR = output.V;
    model.EZ = output.S;
    model.EZZt = zeros(M, M, N);
    for n = 1 : N
        model.EZZt(:,:,n) = output.Sv{n};
    end
else
%  PX  : Noise precision parameters (data)   / point estimate on hyperparam
%  PW  : Noise precision parameters (weight) / point estimate on hyperparam
%  PMU : Noise precision parameters (mean)   / point estimate on hyperparam
%  mZ  : M x N matrix, mean of N latent vectors
%  vZ  : M x M x N cube, variance of N latent vectors
%  mW  : D x M matrix, mean of W
%  vW  : D x M matrix, variance of W
%  mMU : D x 1 vector, mean of Mu
%  vMU : D x 1 matrix, variance of Mu
    model.PX = output.V;
    model.PW = output.Va(1); % WE ASSUMED ALL PRECISIONS TO BE THE SAME
    model.PMU = output.Vmu;
    model.mZ = output.S;
    model.vZ = zeros(M, M, N);
    for n = 1 : N
        model.vZ(:,:,n) = output.Sv{n};
    end
    model.mMU = output.Mu;
    model.vMU = output.Muv;
    model.mW = output.A;
    model.vW = zeros(D, M);
    for d = 1 : D
        for m = 1 : M
            model.vW(d,m) = output.Av{d}(m,m);
        end
    end
end

% cost function, number of iterations, elapsed time
%  eITER    : Iterations took
%  eTIME    : Elapsed time
%  objArray : Objective function value change over iterations
model.objArray = output.lc.cost(2:end);
model.eITER = length(model.objArray);
model.eTIME = output.lc.time(end);

end
