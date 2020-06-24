function model = parfor_cppca( desc_model, varargin  )

fn_model = varargin{1};
X        = varargin{2};
if isequal(desc_model, 'affine')
    [~, W, ~, ~, ~, ~] = fn_model(X');
    model.W = W';   
else
    dim_Z    = varargin{3};
    if isequal(desc_model, 'svd')
        model = feval(fn_model, X, dim_Z);
    elseif isequal(desc_model, 'em')
        objprec_expr       = varargin{4};        
        objfreq_expr       = varargin{5};
        init_model         = varargin{6};
        
        model = feval(fn_model, X, dim_Z,               ...
                     'Threshold',  objprec_expr,        ...
                     'ShowObjPer', objfreq_expr,        ...
                     'InitModel',  init_model               );
    else
        error(['[Error] Undefined Type - ' type_model]);
    end    
end

end

