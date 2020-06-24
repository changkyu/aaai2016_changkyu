function [cm5, cm6] = expr_run_vbpca(X, M)

% PPCA (Ilin)
% disp('PPCA (Ilin)');
% opts = struct( ...
%     'algorithm', 'ppca',... % algorithm
%     'maxiters', 1000,... % maximum iterations
%     'cfstop', [ 1 0 1e-5 ], ... % stop if cost not change over some iters.
%     ... % [#iterations_to_monitor, absolute_error, relative_error]
%     'minangle', 0, ... % keep going regardless of subspace angle change
%     'uniquesv', 0 ... % compute only unique covariance matrices of S(:,j)
% );
% output = pca_full( X, M, opts );
% cm5 = output2model(output);

% VBPCA (Ilin)
% disp('VBPCA (Ilin)');
% opts.algorithm = 'vb';
% output = pca_full( X, M, opts );
% cm6 = output2model(output);

cm5 = [];
cm6 = [];

end
