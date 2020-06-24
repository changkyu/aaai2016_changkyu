function MissIDX = sfm_gen_missing_mnar(X, ratio)
% GEN_SFM_MISSING_MNAR  Missing Not At Random matrix for SfM input matrix
%
% Description
%  Generate a band matrix. When the object is rotating, we can reasonably
%  assume that missing value index matrix may look like a band matrix.
%  Other assumptions are possible but our purpose of the demonstrate is to
%  show that our model is resistent to correlated missing values.
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

k = ceil(sqrt(D*N*ratio));

if k >= D-1
    k = D-2;
end

lower = tril(ones(D,N),-(D-k));
higher = triu(ones(D,N),N-k);

MissIDX = lower + higher;

% make SfM consistent...
idx_start_c = find(MissIDX(1,:) == 1);
idr = 1;
for idx = idx_start_c(1):N
    if mod(idr,2) == 1
        MissIDX(idr, idx) = 0;
    end
    idr = idr + 1;
end

idx_start_r = find(MissIDX(:,1) == 1);
idc = 1;
for idx = idx_start_r(1):D
    if mod(idx,2) == 0
        MissIDX(idx, idc) = 0;
    end
    idc = idc + 1;
end

MissIDX(MissIDX == 1) = -1;
MissIDX(MissIDX == 0) = 1;
MissIDX(MissIDX == -1) = 0;

end