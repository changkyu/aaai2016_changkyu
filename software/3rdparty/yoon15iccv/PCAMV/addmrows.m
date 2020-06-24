%  ADDMROWS - Add unobserved rows to PCA solution
%
%  [ A_new, Av_new ] = ADDMROWS( A, Av, Ir, n1x, Va ) adds the columns
%  removed by RMEMPTY to the found parameter A.
%
%  Va is the variance of the prior for A. Va=inf corresponds to
%  no prior for A.
%
%  See also RMEMPTY, ADDMCOLS

%  This software is provided "as is", without warranty of any kind.
%  Alexander Ilin, Tapani Raiko

function [ A_new, Av_new ] = addmrows( A, Av, Ir, n1x, Va );

ncomp = size(A,2);
if nargin < 5
    Va = inf(1,ncomp);
end

% Might not be very efficient
Ir2 = setdiff( 1:n1x, Ir );

if any(isinf(Va))
    A_new = repmat( NaN, n1x, ncomp );
else
    A_new = zeros(n1x, ncomp);
end
A_new(Ir,:) = A;

if nargout >= 2 & ~isempty(Av)
    if iscell(Av)
        Av_new = cell(1,n1x);
        cur = 1;
        for i = 1:length(Ir)
            Av_new{Ir(i)} = Av{i};
        end
        for i = 1:length(Ir2)
            Av_new{Ir2(i)} = diag(Va);
        end
    else
        % Av is a matrix, each column contains the diagonal
        % elements of the covariance matrix of the corresponding
        % column of A
        
        % Va must be a row vector, each element is a prior for
        % the corresponding column of A
        
        Av_new = repmat( Va, n1x, 1 );
        Av_new(Ir,:) = Av;
    end
else
    Av_new = [];
end
