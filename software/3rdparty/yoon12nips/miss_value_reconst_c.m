function [Xbar] = miss_value_reconst_c(W, EZ, MU)
% MISS_VALUE_RECONSTRUCT_C  Reconstruct input from latent space (central)
%
% INPUT
%   W          Projection matrix
%   EZ         Expected values of latent space
%   MU         Estimated mean
%
% OUTPUT
%   Xbar       Reconstructred input matrix
%
% Implemented/Modified
%  by     Sejong Yoon (sjyoon@cs.rutgers.edu)
%  on     2012.06.01 (last modified on 2012/06/01)

D = size(W, 1);
N = size(EZ, 2);

Xbar = zeros(D, N);

for idd = 1:D
    for idn = 1:N
        Xbar(idd,idn) = W(idd,:) * EZ(:,idn) + MU(idd);
    end
end

end
