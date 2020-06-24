%  ADDMCOLS - Add unobserved columns to PCA solution
%
%  [ S_new, Sv_new, Isv_new ] = ADDMCOLS( S, Sv, Ic, n2x, Isv ) adds
%  the columns removed by RMEMPTY to the found parameter S.
%
%  See also RMEMPTY, ADDMROWS

%  This software is provided "as is", without warranty of any kind.
%  Alexander Ilin, Tapani Raiko

function [ S_new, Sv_new, Isv_new ] = addmcols( S, Sv, Ic, n2x, Isv );

ncomp = size(S,1);
% Might not be very efficient
Ic2 = setdiff( 1:n2x, Ic );

S_new = zeros(ncomp, n2x);
S_new(:,Ic) = S;

if iscell(Sv)
    if nargin < 5 | isempty(Isv)
        Sv_new = cell(1,n2x);
        for j = 1:length(Ic)
            Sv_new{Ic(j)} = Sv{j};
        end
        for j = 1:length(Ic2)
            Sv_new{Ic2(j)} = eye(ncomp);
        end
        Isv_new = [];
    else
        Sv_new = Sv;
        Sv_new{length(Sv)+1} = eye(ncomp);
        Isv_new = repmat( length(Sv_new), 1, n2x );
        Isv_new(Ic) = Isv;
    end
else
    % Sv is a vector of the diagonal elements of a covariance
    % matrix
    Sv_new = repmat( 1, size(S,1), n2x );
    Sv_new(:,Ic) = Sv;
end
    
