function [Xbar] = miss_value_reconst_d(W, EZ, MU)
% MISS_VALUE_RECONSTRUCT_C  Reconstruct input from latent space (dist.)
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

J = size(W, 3);

Nl = zeros(J,1);
Nh = zeros(J,1);
for idj = 1:J
    if idj == 1
        Nl(idj) = 1;
        Nh(idj) = size(EZ{idj},2);
    else
        Nl(idj) = Nh(idj-1) + 1;
        Nh(idj) = Nh(idj-1) + size(EZ{idj}, 2);
    end
end

D = size(W, 1);
N = Nh(J);

Xbar = zeros(D, N);

for idj = 1:J
    Wc = W(:,:,idj);
    MUc = MU(:,idj);
    EZc = EZ{idj};
    Xbarc = zeros(D, Nh(idj)-Nl(idj)+1);

    for idd = 1:D
        for idn = 1:(Nh(idj) - Nl(idj) + 1)
            Xbarc(idd,idn) = Wc(idd,:) * EZc(:,idn) + MUc(idd);
        end
    end
    
    Xbar(:,Nl(idj):Nh(idj)) = Xbarc;
end

end
