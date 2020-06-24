%  MISCOMB - Find combinations of missing values in columns
%
%  [ nobscomb, obscombj, Isv ] = MISCOMB(M,verbose) computes the
%  combinations of missing values in each column of a data matrix.
%  This is needed for faster implementation.
%
%  See also PCA_FULL, PCA_PT, LSPCA

%  This software is provided "as is", without warranty of any kind.
%  Alexander Ilin, Tapani Raiko

function [ nobscomb, obscombj, Isv ] = miscomb(M,verbose)

n2 = size(M,2);

if verbose
    fprintf( 'Calculating Isv ... ');
end
[tmp, I, Isv] = unique( M', 'rows' ); % Transpose of M may be heavy
nobscomb = size(tmp,1); % or length(I) % or max(Isv);

if nobscomb < n2
    obscombj = cell(nobscomb, 1);
    for i = 1:n2
        obscombj{Isv(i)} = [obscombj{Isv(i)}; i];
    end
else
    obscombj = {};
    Isv = [];
end
if verbose
    fprintf('done\n');
    fprintf( 'Missing values combinations: found %d in %d columns\n',...
             nobscomb, n2 )
end
