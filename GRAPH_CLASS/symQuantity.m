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
        function obj = symQuantity(sym, sym_params, args)
            arguments
                sym
                sym_params
                args cell = {}
            end
            
            obj.Sym = sym;
            obj.SymParams = sym_params;
            obj.mtlbFunc = matlabFunction(sym_params,sym,args);
        end
        
        function val = double(obj, sym_param_vals, varargin)
            if nargin == 1
                val = obj.mtlbFunc(obj.SymParams.Vals);
            elseif nargin > 1
                if isempty(sym_param_vals)
                    sym_param_vals = obj.SymParams.Vals;
                else
                    assert(numel(sym_param_vals) == obj.SymParams.N, "Argument requires %d elements", obj.SymParams.N);
                end
                val = obj.mtlbFunc(sym_param_vals,varargin{:});                
            end
        end
    end
end

