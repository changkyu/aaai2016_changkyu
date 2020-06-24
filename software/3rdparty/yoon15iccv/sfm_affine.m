function [sfm_M, sfm_S, sfm_b, U3, W3, V3] = sfm_affine(measurement_matrix)
% SFM_AFFINE  Affine Structure from Motion given measurement matrix. 
%
% Description
%  We first translate points to their center of mass of points in each
% frame, i.e.
%             (x,y) = (x,y) - mean(measurement_matrix,2)
% Then we compute SVD of measurement_matrix as measurement_matrix = U W V'
% Our affine solution for SfM is
%             U   = motion (A)
%             VW' = structure (P)
% We ignore rotational ambiguity for the time being.
%
% IN
%  measurement_matrix - 2M x N measurement matrix with M farmes of N points
%
% OUT
%  sfm_M    - Motion (rotation) matrix (U * sqrt(W) * Q)
%  sfm_S    - Structure matrix (Q' * sqrt(W) * V)
%  sfm_b    - Vector b
%  U3,W3,V3 - Result of SVD on measurement_matrix
%  Q        - Orthogonal constraint

% Compute centroid and set fixed coordinate
N = size(measurement_matrix, 2);

centroid = mean(measurement_matrix, 2);
mm_trans = bsxfun(@minus, measurement_matrix, centroid);

% Do SVD
[U, W, V] = svd(mm_trans);

% We throw away null space
U3 = U(:, 1:3);
W3 = W(1:3, 1:3);
V3 = V(:, 1:3);

% Find an affine SfM solution
sfm_M = U3 * sqrt(W3);
sfm_S = sqrt(W3) * V3';

% Remove ambiguity by finding nullspace (= solving the least square)
[sfm_M,sfm_S] = sfm_remove_rot_amb( sfm_M, sfm_S );

% Reconstruct back
sfm_b = sfm_M * sfm_S + repmat(centroid, [1 N]);

end
