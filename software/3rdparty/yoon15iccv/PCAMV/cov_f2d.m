%  COV_F2D - Extract posterior variance of principal components
%
%  Sv = COV_F2D( Sv, Isv ) extracts the posterior variances of
%  principal components (elements of matrix S) from the structures
%  returned by PCA_FULL.
%
%  See also SUBSPACE2D, ADDEBARS

%  This software is provided "as is", without warranty of any kind.
%  Alexander Ilin, Tapani Raiko

function Sv = cov_f2d( Sv, Isv )

if iscell(Sv)
    Sv_in = Sv;
    if nargin > 1 && ~isempty(Isv)
        Sv = zeros( size(Sv_in{1},1), length(Isv) );
        for i = 1:size(Sv,2)
            Sv(:,i) = diag(Sv_in{Isv(i)});
        end
    else
        Sv = zeros( size(Sv_in{1},1), length(Sv) );
        for i = 1:size(Sv,2)
            Sv(:,i) = diag(Sv_in{i});
        end
    end
end

