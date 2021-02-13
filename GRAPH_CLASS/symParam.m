classdef symParam < sym
    properties
        Default_Value double
    end
    
    methods
        function obj = symParam(sym_arg, def, varargin)            
            if nargin <= 1
                def = 0;
            end
            
            if nargin == 0
                sym_arg = [];
            end
            
            obj@sym(sym_arg, ["real", "positive", varargin{:}]);
            obj.Default_Value = def;

        end
        
        function d = double(obj)
            d = obj.Default_Value;
        end
        
        function d = getDefault(obj)
            d = obj.Default_Value;
        end
    end
end

