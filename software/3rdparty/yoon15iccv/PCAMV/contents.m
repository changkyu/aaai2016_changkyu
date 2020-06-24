% Probabilistic models for principal component analysis (PCA) in the
% presence of missing values. The model for each column of data matrix
% is X(:,j) = Mu + A*S(:,j) + Noise. The noise is isotropic and
% Gaussian with the variance V.
%
% PCA for data sets with missing values:
%   pca_full   - Unrestricted Gaussian posterior for A(i,:) and S(:,j)
%   pca_diag   - Fully factorial approximation of the posterior
%   pca_pt     - Maximum a posteriori estimation
%   pcaimput   - Imputation algorithm
%   lspca      - Least-squares approach
%
% Cost function calculations:
%   cf_full         - Cost for model with unrestricted posterior
%   cf_diag         - Cost for model with fully factorial posterior
%   cf_pt           - Cost for model with MAP estimation
%   compute_rms     - RMS error and reconstruction error matrix
%   errpca_pt.cpp   - Sparse matrix of reconstruction errors
%   errpca_diag.cpp - Sparse matrix of reconstruction errors and
%                     extra parameters needed by PCA_DIAG
%
% Helpers:
%   rmempty         - Remove empty columns or rows from data matrix
%   addmcols        - Add unobserved columns to PCA solution
%   addmrows        - Add unobserved rows to PCA solution
%   miscomb         - Find combinations of missing values in columns
%   argschk         - Check of parameter/value pairs
%   subtract_mu.cpp - Subtract bias term from sparse data matrix
%   converg_check   - Check convergence criteria
%   pcaresults      - Load results from an auto-save file
%
% Visualization:
%   subspace2d - Visualize principal components in 2D
%   covdg      - Extract posterior variance of principal components
%   tsplot     - Plot time series
%   addtsplot  - Add time series to current figure
%   addebars   - Add error bars to current figure
%

%  This software is provided "as is", without warranty of any kind.
%  Alexander Ilin, Tapani Raiko
