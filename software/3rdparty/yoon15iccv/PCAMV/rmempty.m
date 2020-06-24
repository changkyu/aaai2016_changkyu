%  RMEMPTY - Remove empty columns or rows from data matrix
%
%  [ X, Xprobe, Ir, Ic, init ] = RMEMPTY( X, Xprobe, init, verbose )
%  removes empty columns or rows from data matrix X and matrix Xprobe
%  of probing data. Those columns and rows do not affect the found
%  solution.
%
%  See also ADDMCOLS, ADDMROWS

%  This software is provided "as is", without warranty of any kind.
%  Alexander Ilin, Tapani Raiko

function [ X, Xprobe, Ir, Ic, init ] = rmempty( X, Xprobe, init, verbose )

if verbose == 2
    fprintf( 'Checking for empty rows or columns ...\n   ' );
end
    
[n1x,n2x] = size(X);
if issparse(X)
    Ic = find( sum( X~=0, 1 ) > 0 );
    Ir = find( sum( X~=0, 2 ) > 0 );
else
    Ic = find( sum( ~isnan(X), 1 ) > 0 );
    Ir = find( sum( ~isnan(X), 2 ) > 0 );
end
n1 = length(Ir);
n2 = length(Ic);

if n1 == n1x & n2 == n2x
    if verbose
        fprintf( 'No empty rows or columns\n' );
    end
    return
end

%X = X(Ir,Ic);
%if ~isempty(Xprobe), Xprobe = Xprobe(Ir,Ic); end

if n1 < n1x && n2 < n2x
    X = X(Ir,Ic);
    if ~isempty(Xprobe), Xprobe = Xprobe(Ir,Ic); end

elseif n1 < n1x && n2 == n2x
    X = X(Ir,:);
    if ~isempty(Xprobe), Xprobe = Xprobe(Ir,:); end

elseif n1 == n1x && n2 < n2x
    X = X(:,Ic);
    if ~isempty(Xprobe), Xprobe = Xprobe(:,Ic); end

else
    Ir = []; Ic = [];
end

if verbose
    fprintf( '%d empty rows and %d empty columns removed\n',...
             n1x-n1, n2x-n2 );
end

if nargout < 5 || ~isstruct(init)
    return
end

if n1 < n1x
    if isfield( init, 'A' ), init.A = init.A(Ir,:); end
    
    if isfield( init, 'Av' ) & ~isempty(init.Av)
        if iscell( init.Av )
            Av = cell(1,n1);
            for i = 1:n1
                Av{i} = init.Av{Ir(i)};
            end
            init.Av = Av;
            clear Av
        else
            init.Av = init.Av(Ir,:);
        end
    end
    
    if isfield( init, 'Mu' ), init.Mu = init.Mu(Ir,:); end
    
    if isfield( init, 'Muv' ) & ~isempty(init.Muv)
        init.Muv = init.Muv(Ir,:);
    end

end

if n2 < n2x
    if isfield( init, 'S' ), init.S = init.S(:,Ic); end
    
    if isfield( init, 'Sv' ) & ~isempty(init.Sv)
        if iscell( init.Sv )
            if ~isfield( init, 'Isv' ) || isempty(init.Isv)
                Sv = cell(1,n2);
                for j = 1:n2
                    Sv{j} = init.Sv{Ic(j)};
                end
                init.Sv = Sv;
            else
                init.Isv = init.Isv(Ic);
                [B,I,J] = unique( init.Isv );
                init.Sv = {init.Sv{B}};
                init.Isv = J;
            end
        else
            init.Sv = init.Sv(:,Ic);
        end
    end
end
