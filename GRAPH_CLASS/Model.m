classdef Model < matlab.mixin.Copyable
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
        LinearizeFlag logical = true % Flag that determines whether to calculate the linear models
    end

    properties (SetAccess = protected) 
        Nx % number of states
        Nu % number of inputs
        Nd % number of disturbances
        Ny % number of outputs
        f_sym (:,1) sym % f(x,u,d), can contain symbolic parameters
        g_sym (:,1) sym % g(x,u,d), can contain symbolic parameters

        SymVars struct = struct() % Contains fields x,u,d, each an array of symbolic variables

        SymParams sym = sym.empty()
        SymParams_Vals double = []
        N_SymParams double = [] % It takes a bunch of time to count symbolic variables
        
        LinModel LinearModel = LinearModel.empty() % this could be another object
    end
    
    properties (SetAccess = protected, GetAccess = protected)
        CalcF_func function_handle {mustBeScalarOrEmpty} % calculates x_dot
        CalcG_func function_handle {mustBeScalarOrEmpty} % calculates y   
    end
    
    properties (Dependent)
        StateNames (:,1) string
        InputNames (:,1) string
        DisturbanceNames (:,1) string
        OutputNames (:,1) string    
    end
     
    methods
        function obj = Model(varargin)
            % add functionality here at some point
%             if nargin == 0
%                 % do nothing
%             elseif nargin == 1
%                 disp(sprintf('Model of class %s created.\n',class(varargin{1}))) % update this to list valid arguments
% 
%             elseif nargin == 2
%                 obj.f_sym = varargin{1};
%                 obj.g_sym = varargin{2};
%                 init(obj);
%             end    
        end
        
        function init(obj)
            setSymVars(obj);
            SymbolicSolve(obj);
            setCalcFuncs(obj);

            if obj.LinearizeFlag
                obj.setLinearModel;
            end
            
            obj.N_SymParams = numel(obj.SymParams);
        end

        function setSymVars(obj)
            x1 = genSymVars('x%d', obj.Nx); % Symbolic Vars for Dynamic States
            u1 = genSymVars('u%d', obj.Nu); % Symbolic Vars for Inputs
            d1 = genSymVars('d%d', obj.Nd); % Symbolic Vars for Disturbances
            
            obj.SymVars.x = x1;
            obj.SymVars.u = u1;
            obj.SymVars.d = d1;           
        end
        
        function setCalcFuncs(obj)
            f = obj.f_sym;
            g = obj.g_sym;
            
            vars = {[obj.SymVars.x], [obj.SymVars.u], [obj.SymVars.d]};
            CalcFuncs_Cell = genMatlabFunctions(obj, {f, g},vars);
            obj.CalcF_func = CalcFuncs_Cell{1};
            obj.CalcG_func = CalcFuncs_Cell{2};
        end
        
        function setLinearModel(obj)   
            f = obj.f_sym;
            g = obj.g_sym;
            
            A = jacobian(f,obj.SymVars.x);
            B = jacobian(f,obj.SymVars.u);
            E = jacobian(f,obj.SymVars.d);
            
            C = jacobian(g,obj.SymVars.x);
            D = jacobian(g,obj.SymVars.u);
            H = jacobian(g,obj.SymVars.d);
            
            obj.LinModel = LinearModel(A,B,E,C,D,H);
            
            % obj.LinModel.CalcState  = matlabFunction(obj.LinModel.A_sym,obj.LinModel.B_sym,obj.LinModel.E_sym,'Vars',vars); - This should go in LinModel
            % obj.LinModel.CalcOutput = matlabFunction(obj.LinModel.C_sym,obj.LinModel.D_sym,obj.LinModel.H_sym,'Vars',vars);
            %             obj.LinearModel.CalcOutput = matlabFunction(obj.LinearModel.B_sym,'Vars',[{[x1] [u1], [d1]}]);
            %             obj.LinearModel.CalcE = matlabFunction(obj.LinearModel.E_sym,'Vars',[{[x1] [u1], [d1]}]);
            %             obj.LinearModel.CalcC = matlabFunction(obj.LinearModel.C_sym,'Vars',[{[x1] [u1], [d1]}]);
            %             obj.LinearModel.CalcD = matlabFunction(obj.LinearModel.D_sym,'Vars',[{[x1] [u1], [d1]}]);
            %             obj.LinearModel.CalcH = matlabFunction(obj.LinearModel.H_sym,'Vars',[{[x1] [u1], [d1]}]);
        end
        
        function F = CalcF(obj,x,u,d,params)
            param_lengths = [obj.Nx, obj.Nu, obj.Nd, obj.N_SymParams];
            
            if nargin == 4 || obj.N_SymParams == 0
                vars = {x,u,d};
            elseif nargin == 5
                vars = {x,u,d,params};
            end
            
            for i = 1:numel(vars)
                assert(size(vars{i},1) >= param_lengths(i), "Argument %d requires %d entries", i, param_lengths(i));
            end
            
            F = CalcX(obj, obj.CalcF_func, vars);
        end
         
        function G = CalcG(obj,x,u,d,params)
            param_lengths = [obj.Nx, obj.Nu, obj.Nd, obj.N_SymParams];
            
            if nargin == 4 || obj.N_SymParams == 0
                vars = {x,u,d};
            elseif nargin == 5
                vars = {x,u,d,params};
            end
            
            for i = 1:numel(vars)
                assert(size(vars{i},1) >= param_lengths(i), "Argument %d requires %d entries", i, param_lengths(i));
            end
            
            G = CalcX(obj, obj.CalcG_func, vars);
        end
        
        function [A,B,E,F0,C,D,H,G0] = Linearize(obj,x0,u0,d0)
            [A,B,E] = obj.LinModel.CalcState(x0,u0,d0);
            F0 = obj.CalcF_func(x0,u0,d0) - A*x0 - B*u0 - E*d0;
            
            [C,D,G]  = obj.LinModel.CalcOutput(x0,u0,d0);
            G0 = obj.CalcG_func(x0,u0,d0) - C*x0 - D*u0 - G*d0;     
        end
        
        function x = get.StateNames(obj)
            x = defineStateNames(obj);
        end
        
        function x = get.InputNames(obj)
            x = defineInputNames(obj);
        end
        
        function x = get.DisturbanceNames(obj)
            x = defineDisturbanceNames(obj);                       
        end
        
        function x = get.OutputNames(obj)
            x = defineOutputNames(obj);                       
        end     
    end
    
    methods (Access = protected)
        function funcs = genMatlabFunctions(obj, syms, vars)
            % Generates matlabFunctions from symbolic arrays
            % i.e. f_sym -> calcF_Func and g_sym -> calcG_Func in 
            % setCalcFuncs().
            
            if ~isempty(obj.SymParams)
                vars{end+1} = [obj.SymParams];
            end
            
            cell_flag = isa(syms, 'cell');
            if cell_flag
                funcs = cellfun(@(x) matlabFunction(x,'Vars',vars), syms, 'UniformOutput', false);
            else
                funcs = matlabFunction(syms,'Vars',vars);
            end
        end
        
        function X = CalcX(obj, func, vars)
            % Wrapper for matlabFunction properties, does error checking 
            % and assists with vectorizing the function
            
            n_ins = nargin(func);
            assert(numel(vars) == n_ins, "Func Requires %d Arguments", n_ins);
                        
            if size(vars{1},2) == 1
                X = func(vars{:});
            else
                X = splitapply(func,vars{:},1:size(vars{1},2));
            end
        end
    end
    
    methods (Abstract)
        % The following abstract methods must be defined in 'Model's' subclasses
        
        SymbolicSolve % Required to calculate f_sym and g_sym
  
        defineStateNames
        defineInputNames
        defineDisturbanceNames
        defineOutputNames
    end
end

