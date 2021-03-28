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
    % Add Constructor
    % Find better solution than splitapply() for calcX
    % VPA all the things
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    properties
        Nx (1,1) double = 0 % number of states
        Nu (1,1) double = 0 % number of inputs
        Nd (1,1) double = 0 % number of disturbances
        Ny (1,1) double = 0 % number of outputs
        
        StateDescriptions string
        InputDescriptions string
        DisturbanceDescriptions string
        OutputDescriptions string
        
        SymVars SymVars {mustBeScalarOrEmpty} % Contains fields x,u,d, each an array of symbolic variables
        SymParams SymParams {mustBeScalarOrEmpty}
    end
    
    properties
        f_sym (:,1) sym % f(x,u,d), can contain symbolic parameters
        g_sym (:,1) sym % g(x,u,d), can contain symbolic parameters
    end
    
    properties (SetAccess = protected, GetAccess = protected)
        f_func function_handle {mustBeScalarOrEmpty} % calculates x_dot
        g_func function_handle {mustBeScalarOrEmpty} % calculates y
    end
    
    properties (Dependent)
        StateTable table
        InputTable table
        DisturbanceTable table
        OutputTable table
    end
    
    methods
        function init(obj)
            if isempty(obj.SymVars)
                setSymVars(obj);
            end
            setCalcFuncs(obj);
        end
        
        function setSymVars(obj)
            obj.SymVars = SymVars('Nx', obj.Nx, 'Nd', obj.Nd, 'Nu', obj.Nu);
        end
        
        function setCalcFuncs(obj)
            f = obj.f_sym;
            g = obj.g_sym;
            
            CalcFuncs_Cell = genMatlabFunctions(obj, {f, g});
            obj.f_func = CalcFuncs_Cell{1};
            obj.g_func = CalcFuncs_Cell{2};
        end
        
        function F = CalcF(obj,x,u,d,params)
            param_lengths = [obj.Nx, obj.Nu, obj.Nd, obj.SymParams.N];
            
            if isempty(obj.SymParams)
                vars = {x,u,d};
            else
                if nargin == 4
                    vars = {x,u,d,obj.SymParams.Vals};
                elseif nargin == 5
                    vars = {x,u,d,params};
                end
            end
            
            for i = 1:numel(vars)
                assert(size(vars{i},1) >= param_lengths(i), "Argument %d requires %d entries", i, param_lengths(i));
            end
            
            F = obj.CalcX(obj.f_func, vars);
        end
        
        function G = CalcG(obj,x,u,d,params)
            param_lengths = [obj.Nx, obj.Nu, obj.Nd, obj.SymParams.N];
            
            if isempty(obj.SymParams)
                vars = {x,u,d};
            else
                if nargin == 4
                    vars = {x,u,d,obj.SymParams.Vals};
                elseif nargin == 5
                    vars = {x,u,d,params};
                end
            end
            
            for i = 1:numel(vars)
                assert(size(vars{i},1) >= param_lengths(i), "Argument %d requires %d entries", i, param_lengths(i));
            end
            
            G = obj.CalcX(obj.g_func, vars);
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
        
        function t = get.StateTable(obj)
            state_syms = arrayfun(@(x) string(sym2str(x)), obj.SymVars.x);
            t = table(state_syms, obj.StateDescriptions);
        end
        
        function t = get.InputTable(obj)
            in_syms = arrayfun(@(x) string(sym2str(x)), obj.SymVars.u);
            t = table(in_syms, obj.InputDescriptions);
        end
        
        function t = get.DisturbanceTable(obj)
            dist_syms = arrayfun(@(x) string(sym2str(x)), obj.SymVars.d);
            t = table(dist_syms, obj.DisturbanceDescriptions);
        end
        
        function t = get.OutputTable(obj)
            out_syms = arrayfun(@(x) sprintf("y%d", x), 1:obj.Ny);
            t = table(out_syms', obj.OutputDescriptions);
        end
    end
   
    methods (Access = protected)
        function funcs = genMatlabFunctions(obj, syms, vars)
            % Generates matlabFunctions from symbolic arrays
            % i.e. f_sym -> calcF_Func and g_sym -> calcG_Func in
            % setCalcFuncs().
            % Vars argument is optional
            arguments 
                obj
                syms
                vars = {}
            end
            
            if isempty(vars)
                vars = {[obj.SymVars.x], [obj.SymVars.u], [obj.SymVars.d]};
            end
            
            if ~isempty(obj.SymParams)
                vars{end+1} = [obj.SymParams.Syms];
            end
            
            cell_flag = isa(syms, 'cell');
            if cell_flag
                funcs = cellfun(@processSym, syms, 'UniformOutput', false);
            else
                funcs = processSym(syms);
            end
            
            function func = processSym(sym)
                if isa(sym, 'sym')
                    func = matlabFunction(sym,'Vars',vars);
                else
                    func = @(varargin) sym;
                end
            end
        end
        
        function copyModelProps(obj_from, obj_to, opts)
            arguments
                obj_from
                obj_to
                opts.Properties = ["Nx","Nu","Nd","Ny",...
                "StateDescriptions", "InputDescriptions", "DisturbanceDescriptions", "OutputDescriptions",...
                "SymVars","SymParams"];
            end

            for prop = opts.Properties
                obj_to.(prop) = obj_from.(prop);
            end
        end
    end
    
    methods (Static, Access = protected)
        function X = CalcX(func, vars)
            % Wrapper for matlabFunction properties, does error checking
            % and assists with vectorizing the function
            
            n_ins = nargin(func);
            if n_ins > -1
                assert(numel(vars) == n_ins, "Func Requires %d Arguments", n_ins);
            end
            
            if size(vars{1},2) == 1
                X = func(vars{:});
            else
                X = splitapply(func,vars{:},1:size(vars{1},2));
            end
        end
    end
end

