function model = cppca( X, M )
% CPPCA      Probablistic Principal Component Analysis (PPCA)
% 
% Description
%  model = ppca(X, M) computes MLE closed form solution of W and VAR given 
%  (D dimension)x(N observations) sample matrix X. Technically identical 
%  but numerically slightly different from PMTK3 (ppcaFit.m). This is a
%  reference implementation and does not handle missing values.
%
% Output
%  model = structure(W, MU, VAR, U, EZ, EZZt, Lambda)
%  W     : D x M projection matrix
%  MU    : D x 1 vector sample means
%  VAR   : Scalar estimated variance
%  U     : D x M matrix containing column vectors of eigenvectors
%  EZ    : M x N matrix latent space mean
%  EZZt  : M x M matrix latent space covariance matrix
%  Lambda: D x 1 vector singular values
%
% Implemented
%  by     Sejong Yoon (sjyoon@cs.rutgers.edu)
%  on     2011.10.05
%
% References
%  [1] M.E. Tipping and C.M. Bishop, Probablistic principal component 
%      analysis, J. R. Statist. Soc. B (1999), 61 Part 3, pp. 611-622.
%  [2] Probablistic Modeling Toolkit 3, pmtk3.googlecode.com

% D = original dimension
[D, N] = size(X);

% Compute sample mean and covariance matrix
MU = mean(X, 2);

% Compute eigenvectors and eigenvalues of covariance matrix
% Compute only needed singular values since M < D in many cases
[~, S, U] = svd((X - repmat(MU, [1 N]))',0);
Lambda = diag(S).^2;

% Find principal components
[Lambda, idx] = sort(Lambda, 'descend');
U = U(:,idx);

% VAR=(1/D-M)*(LAMBDA(M+1)+...+LAMBDA(D))
VAR = sum(Lambda(M+1:end)) / (D - M);

% W=U*{(LAMBDA - VAR*I)^(1/2)}*R
W = U(:,1:M) * (diag(Lambda(1:M)) - VAR * eye(M))^(0.5);

% Compute E[Z] and E[ZZ']
% E[Z]   = M^-1 * W' * (X - MU)
% E[ZZ'] = sigma^2 * M^-1 + E[Z]E[Z]'
%   where M = W'W + sigma^2*I
EZ = (W'*W + VAR*eye(M))^(-1) * W' * (X - repmat(MU, [1 N]));
EZZt = VAR * (W'*W + VAR*eye(M))^(-1) + EZ * EZ';

% Create structure
model = structure(W, MU, VAR, U, EZ, EZZt, Lambda);

end
