classdef symParam2 < handle
    properties
        Sym sym = sym.empty()
        Default_Value double = []
    end
    
    properties (SetAccess = ?Component)
        Parent Component = Component.empty()
    end
    
    methods
        function obj = symParam2(sym_arg, def, varargin)            
            if nargin <= 1
                def = 0;
            end
            
            if nargin == 0
                sym_arg = [];
            end
            
            obj.Sym = sym(sym_arg, varargin{:});
            obj.Default_Value = def;
        end
        
        function d = double(obj)
            d = obj.Default_Value;
        end
        
        function s = sym(obj)
            s = obj.Sym;
        end
    end
end

