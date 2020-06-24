function [model_init] = get_init_value_ex(X, dim_Z, varargin)
    if( numel(varargin)== 1 )
        model_init = get_init_value(X, dim_Z, varargin{1}{:});
    else
        model_init = get_init_value(X, dim_Z, varargin{:});
    end
end

