function [MissIDX] = gen_sfm_missing_mar(X, ratio)
% GEN_SFM_MISSING_MAR   Missing At Random matrix for SfM input matrix
%
% Description
%  Generate missing value index matrix with missing-at-random
%
% Input
% X        : Input matrix
% ratio    : Missing ratio
%
% Output
% MissIDX  : Missing value matrix with same size as X
%
% Implemented
%  by     Sejong Yoon (sjyoon@cs.rutgers.edu)
%  on     2011.10.07 (last modified on 2012/02/01)

[D, N] = size(X);

if mod(D,2) ~= 0
    error('Features should be twice the number of frames!');
elseif ratio >= 1 || ratio < 0
    error('Ratio should be in range [0 1)!');
else
    F = D / 2;
end

while 1
    % we randomly miss some points
    MissIDX = ones(1, N*F);
    MissIDX(randperm(N*F, ceil(ratio * N*F))) = 0;
    MissIDX = reshape(MissIDX, [N F]);
    MissIDX = repmat(MissIDX, [2 1]);
    MissIDX = reshape(MissIDX, [N F*2])';

    % check if at least one frame is observed for each point
    if sum(sum(MissIDX, 1) == 0) == 0
        break;
    end
end

% flip
MissIDX = ~MissIDX;

end
