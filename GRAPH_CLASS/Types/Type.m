classdef Type < matlab.mixin.Copyable
    %TYPE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties % User can set Val_Char or Val_Sym, the corresponding properties are updated automatically
        Val_Char string = string.empty()
        Val_Sym sym = sym.empty()
    end
    
    properties (SetAccess = protected)
        vars sym = sym.empty % contains list of symbolic variables used in type definition.  This property is set by subclasses, i.e. Type_PowerFlow.vars = [x_t, x_h, u_j]
        
        Val_Func function_handle = function_handle.empty() % Value calculation function 
        Jac_Func function_handle = function_handle.empty() % Jacobian calculation function
        
        Jac_Sym sym = sym.empty() % Symbolic Jacobian
    end
    
    methods
        function obj = Type(varargin)
            if nargin > 1
                obj = my_inputparser(obj,varargin{:});
            end
        end
        
        function set.Val_Char(obj, val)
            obj.Val_Char = val;
        end
        
        function set.Val_Sym(obj, val)
        end
        
        function val = calcVal(obj, vars_) % Calculates type value with symbolic 'vars' substituted with numeric 'vars_'
        end
        
        function jac = calcJac(obj, vars_) % Calculates type Jacobian with symbolic 'vars' substituted with numeric 'vars_'
        end
        
    end
end

