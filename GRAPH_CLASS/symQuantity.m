classdef symQuantity
    % Stores symbolic expressions composed of SymParams elements
 
    properties
        Sym {mustBeNumericOrSym}
        SymParams SymParams
    end
    
    properties (SetAccess = private, GetAccess = private)
        mtlbFunc function_handle
    end
    
    methods
        function obj = symQuantity(sym, sym_params)
            obj.Sym = sym;
            obj.SymParams = sym_params;
            obj.mtlbFunc = matlabFunction(sym_params,sym);
        end
        
        function val = double(obj, sym_param_vals)
            if nargin == 1
                val = obj.mtlbFunc(obj.SymParams.Vals);
            elseif nargin == 2
                assert(numel(sym_param_vals) == obj.SymParams.N, "Argument requires %d elements", obj.SymParams.N);
                val = obj.mtlbFunc(sym_param_vals);
            end
        end
    end
end

