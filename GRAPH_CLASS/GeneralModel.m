classdef GeneralModel < Model_Super
    % The Model class in the Graph Modeling Toolbox is used to generically
    % define a model in nonlinear state space form.
    %
    % System Description
    % x_dot = f_sym(x,u,d)
    % y     = g_sym(x,u,d)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Contributors: Christopher T. Aksland and Phil Renkert
    % Association: University of Illionis at Urbana-Champaign
    % Contact: aksland2@illinois.edu and renkert2@illinois.edu
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % potential improvements:
    % Phil: Figure out more elegant SymParam...
    % Add Constructor
    % Find better solution than splitapply() for calcX
    % VPA all the things
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    properties
        f_sym (:,1) sym % f(x,u,d), can contain symbolic parameters
        g_sym (:,1) sym % g(x,u,d), can contain symbolic parameters
    end
    
    properties (SetAccess = protected, GetAccess = protected)
        f_func function_handle {mustBeScalarOrEmpty} % calculates x_dot
        g_func function_handle {mustBeScalarOrEmpty} % calculates y   
    end
   
    methods
        function init(obj)
            init@Model_Super(obj);
            setCalcFuncs(obj);
        end
        
        function setCalcFuncs(obj)
            f = obj.f_sym;
            g = obj.g_sym;
            
            CalcFuncs_Cell = genMatlabFunctions(obj, {f, g});
            obj.f_func = CalcFuncs_Cell{1};
            obj.g_func = CalcFuncs_Cell{2};
        end
        
        function lm = getLinearModel(obj)   
            f = obj.f_sym;
            g = obj.g_sym;
            
            A = jacobian(f,obj.SymVars.x);
            B = jacobian(f,obj.SymVars.u);
            E = jacobian(f,obj.SymVars.d);
            
            C = jacobian(g,obj.SymVars.x);
            D = jacobian(g,obj.SymVars.u);
            H = jacobian(g,obj.SymVars.d);
            
            lm = LinearModel(A,B,E,C,D,H);
            copyModelProps(obj, lm);
        end
        
        function F = CalcF(obj,x,u,d,params)
            param_lengths = [obj.Nx, obj.Nu, obj.Nd, obj.SymParams.N];
            
            if nargin == 4 || isempty(obj.SymParams)
                vars = {x,u,d};
            elseif nargin == 5
                vars = {x,u,d,params};
            end
            
            for i = 1:numel(vars)
                assert(size(vars{i},1) >= param_lengths(i), "Argument %d requires %d entries", i, param_lengths(i));
            end
            
            F = obj.CalcX(obj.f_func, vars);
        end
         
        function G = CalcG(obj,x,u,d,params)
            param_lengths = [obj.Nx, obj.Nu, obj.Nd, obj.SymParams.N];
            
            if nargin == 4 || isempty(obj.SymParams)
                vars = {x,u,d};
            elseif nargin == 5
                vars = {x,u,d,params};
            end
            
            for i = 1:numel(vars)
                assert(size(vars{i},1) >= param_lengths(i), "Argument %d requires %d entries", i, param_lengths(i));
            end
            
            G = obj.CalcX(obj.g_func, vars);
        end      
    end    
end

