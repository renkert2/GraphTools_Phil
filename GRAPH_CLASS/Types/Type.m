classdef Type < matlab.mixin.Copyable
    % Type is a Super Class for all mathmatical expressions in the the 
    % Graph Modeling Toolbox. For example, objects of class Type can
    % convert an expression, x1+2*x2, into a callable matlab function. Type
    % objects are typically used to define graph powerflows and
    % capacitances, but also have other more generic applications. 
    % Instatiate an empty object or a filled object in the following form:
    % 
    % T = Type(type, params) where 
    % - type is a desired function defined as a string or symbolic expression
    % - params is a cell array of symbolic variables in type used to define the argument structure, 
    % -- When defining the internal matlabFunctions, the params property is used in
    % -- matlabFunction(___, 'Vars', obj.params)
    % EX: T = Type('b1*b2*c', {[sym('b1'); sym('b2')], sym('c')})
    
    % T = Type(type) 
    % When only type is specified, arguments to calcVal are determined by symvar(type)
    % ex: T = Type('b1*b2')
      
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Contributors: Christopher T. Aksland and Phil Renkert
    % Association: University of Illionis at Urbana-Champaign
    % Contact: aksland2@illinois.edu and renkert2@illinois.edu
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % potential improvements:
    % - remove either vars or params since they are duplicates of each
    % other. The redundancy may cause confusion.
    % - add 'optimize' option to the matlabfunction calls to speed up
    % processing?
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    properties  % User can set Val_Str or Val_Sym, the corresponding properties are updated automatically
        Val_Str string = string.empty()
        Val_Sym sym = sym.empty()
    end
    
    properties (SetAccess = protected)
        params cell % Cell array of variables used to define matlabFunction input arguments
        Val_Func function_handle % Value calculation function
    end
    
    properties (SetAccess = private, GetAccess = private)
        set_method_flag logical = true % Allows access to set methods, prevents recursion
    end
    
    methods
        function obj = Type(type, params)
            if nargin ~= 0
                if nargin == 1 % Only type specified
                    setType(obj, type);
                    obj.params = sym2cell(symvar(obj.Val_Sym));
                elseif nargin ==2
                    setType(obj, type);
                    obj.params = params;
                end
                setCalcProps(obj)
            end
        end
        
        function init(obj, type)
            setType(obj, type);
            setCalcProps(obj);
        end
        
        function setType(obj,type)
            if isa(type,'string') || isa(type,'char')
                Val_Str_temp = type;
                Val_Sym_temp = str2sym(Val_Str_temp);
            elseif isa(type,'sym')
                Val_Sym_temp = type;
                Val_Str_temp = string(Val_Sym_temp);
            else
                error("Type must be 'string', 'char', or 'sym'")
            end
                                 
            obj.set_method_flag = false; % Prevents infinite loop in set method
            
            if isempty(obj.Val_Sym) || obj.Val_Sym ~= Val_Sym_temp
                obj.Val_Sym = Val_Sym_temp;
            end
            
            if isempty(obj.Val_Str) || not(strcmp(obj.Val_Str,Val_Str_temp))
                obj.Val_Str = Val_Sym_temp;
            end
            
            obj.set_method_flag = true;
        end
        
        function setCalcProps(obj)
            obj.Val_Func = matlabFunction(obj.Val_Sym,'Vars',obj.params);
        end
        
        % set the string value and update the symbolic definition
        function set.Val_Str(obj, type)
            if obj.set_method_flag
                init(obj, type)
            else
                obj.Val_Str = type;
            end
        end
        
        % set the symbolic value and update the string definition
        function set.Val_Sym(obj, type)
            if obj.set_method_flag
                init(obj,type)
            else
                obj.Val_Sym = type;
            end
        end
        
        
        function val = calcVal(obj, varargin)
            % Calculates type value by replacing symbolic vars with numerical vars.
            % obj can be a single Type object or object array of Types with similar arguments
            % Expects one argument for each element in obj.params.  
            % E.X. if obj.params = {[x1;x2],y,z}, use obj.calcVal([x1;x2], y, z)
            
            num_params = arrayfun(@(x) numel(x.params), obj);
            assert(numel(varargin) == max(num_params), 'Arguments must match Parameter Structure.  For object array, ensure maximum number of inputs are given');
            
            if numel(obj)==1
                val = obj.Val_Func(varargin{1:num_params});
            else
                val = cell(size(obj));
                for i = 1:size(obj,1)
                    for j = 1:size(obj,2)
                        val{i,j} = obj(i,j).Val_Func(varargin{1:num_params(i,j)});
                    end
                end
            end
        end
    end
end
    
