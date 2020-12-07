classdef Type < matlab.mixin.Copyable
    %TYPE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties % User can set one of these two properties, the corresponding property is updated automatically
        Val_Char char = char.empty
        Val_Sym sym = sym.empty
    end
    
    properties (SetAccess = protected)
        Val_Func function_handle % Value calculation function 
        Jac_Func function_handle % Jacobian calculation function
        
        Jac_Sym sym % Symbolic Jacobian
    end
    
    methods
        function val = calcVal(obj)
        end
        
        function jac = calcJac(obj)
        end
        
    end
end

