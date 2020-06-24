function err = calc_ppca_rms( X, W, Z, MU )

if iscell(Z) 
    % distributed case
    Jz = length(Z);
    if ~iscell(X)
        % we assume X and Z are partitioned in the same way
        Xt = cell(Jz, 1);
        ub = 0;
        for i = 1 : Jz
            lb = ub + 1;
            ub = lb + size(Z{i},2) - 1;
            Xt{i} = X(:,lb:ub);
        end
        X = Xt;
    end
    
    Jx = length(X);
    Jw = length(W);
    Jm = length(MU);
    assert(Jx == Jz, 'Invalid Z : dimension mismatch');
    assert(Jw == Jz, 'Invalid W : dimension mismatch');
    assert(Jm == Jz, 'Invalid MU: dimension mismatch');
    
    XWZplusMU = cell(1,Jz);
    for j = 1 : Jz
        Dx = size(X{j}, 1);
        [Dw, Mw] = size(W{j});
        Mz = size(Z{j}, 1);
        Dm = size(MU{j}, 1);
        assert(Mw == Mz, 'Invalid W or Z: dimension mismatch');
        assert(Dx == Dw, 'Invalid X or W: dimension mismatch');
        assert(Dx == Dm, 'Invalid X or MU: dimension mismatch');
        
        XWZplusMU{j} = X{j} - bsxfun(@plus, W{j}*Z{j}, MU{j});
    end
    
    XWZplusMU = cell2mat(XWZplusMU);
    X = cell2mat(X');
else
    % centralized case
    [Dx, Nx] = size(X);
    [Mz, Nz] = size(Z);
    [Dw, Mw] = size(W);
    [Dm, ~] = size(MU);
    assert(Dx == Dw, 'Invalid X or W: dimension mismatch');
    assert(Nx == Nz, 'Invalid X or Z: dimension mismatch');
    assert(Mz == Mw, 'Invalid Z or W: dimension mismatch');
    assert(Dm == Dx, 'Invalid MU or X: dimension mismatch');
    
    XWZplusMU = X - bsxfun(@plus, W*Z, MU);
end

err = rms(rms( XWZplusMU(~isnan(X)) ));

end