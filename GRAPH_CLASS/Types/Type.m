classdef Type < matlab.mixin.Copyable
    %TYPE Summary of this class goes here
    %   Detailed explanation goes here
   
    
    properties % User can set Val_Char or Val_Sym, the corresponding properties are updated automatically
        Val_Char string = string.empty() % Should this be Called Val_Str... Actually, maybe we should name things better for input parsing...?
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
        
        % set the string value and update the symbolic definition
        function set.Val_Char(obj, val) 
            obj.Val_Char = val;
            obj = updateVal(obj);
        end
        
        % set the symbolic value and update the string definition
        function set.Val_Sym(obj, val)
            obj.Val_Sym = val;
            obj = updateVal(obj);
        end
        
        % update the type value property that was not defined
        function obj = updateVal(obj)
            if isempty(obj.Val_Sym)
                obj.Val_Sym = str2sym(obj.Val_Char);
            end
            if isempty(obj.Val_Char)
                obj.Val_Char ='Val_Char', string(obj.Val_Sym);
            end
        end
        
        
        
        function val = calcVal(obj, vars_) % Calculates type value with symbolic 'vars' substituted with numeric 'vars_'
        end
        
        function jac = calcJac(obj, vars_) % Calculates type Jacobian with symbolic 'vars' substituted with numeric 'vars_'
        end
        
    end
end

