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
        
        x0 double
        
        StateDescriptions (:,1) string
        InputDescriptions (:,1) string
        DisturbanceDescriptions (:,1) string
        OutputDescriptions (:,1) string
        
        SymVars SymVars {mustBeScalarOrEmpty} % Contains fields x,u,d, each an array of symbolic variables
        Params compParam
        
        Name string % Used to identify model subsystem blocks in Simulink model
    end
    
    properties
        f_sym (:,1) sym % f(x,u,d), can contain symbolic parameters
        g_sym (:,1) sym % g(x,u,d), can contain symbolic parameters
    end
    
    properties (SetAccess = protected)
        LinearModel LinearModel
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
        
        function set.x0(obj, x0)
            if obj.Nx
                assert(numel(x0) == obj.Nx, 'Value for x0 must have %d entries', obj.Nx);
            end
            obj.x0 = x0;
        end
        
        function F = CalcF(obj,x,u,d)
            param_lengths = [obj.Nx, obj.Nu, obj.Nd];
            
            vars = {x,u,d};
            
            for i = 1:numel(vars)
                assert(size(vars{i},1) >= param_lengths(i), "Argument %d requires %d entries", i, param_lengths(i));
            end
            
            F = obj.CalcX(obj.f_func, vars);
        end
        
        function G = CalcG(obj,x,u,d)
            param_lengths = [obj.Nx, obj.Nu, obj.Nd];
            
            vars = {x,u,d};
            
            for i = 1:numel(vars)
                assert(size(vars{i},1) >= param_lengths(i), "Argument %d requires %d entries", i, param_lengths(i));
            end
            
            G = obj.CalcX(obj.g_func, vars);
        end
        
        function [x_bar,y_bar] = calcSteadyState(obj, u, d, x0, opts)
            % Solves CalcF == 0 to get the steady state values
            arguments
                obj
                u double
                d double
                x0 double = []
                opts.SolverOpts struct = optimset('Display', 'off');
            end
            
            if isempty(x0)
                x0 = zeros(obj.Nx,1);
            end
            [x_bar] = fsolve(@(x) CalcF(obj,x,u,d), x0, opts.SolverOpts);
            
            if nargout == 2
                y_bar = CalcG(obj,x_bar,u,d);
            end
        end
        
        function [u_bar, y_bar] = calcSteadyStateInput(obj, x, d, u0, opts)
            % Solves CalcF == 0 to get the steady state input
            arguments
                obj
                x double
                d double
                u0 double = []
                opts.SolverOpts struct = optimset('Display', 'off');
            end
            
            if isempty(u0)
                u0 = zeros(obj.Nu,1);
            end
            [u_bar] = fsolve(@(u) CalcF(obj,x,u,d), u0, opts.SolverOpts);
            
            if nargout == 2
                y_bar = CalcG(obj,x,u_bar,d);
            end
        end
        
        function lm = getLinearModel(obj)
            f = obj.f_sym;
            g = obj.g_sym;
            
            x = obj.SymVars.x;
            u = obj.SymVars.u;
            d = obj.SymVars.d;
            
            A = jacobian(f,x);
            B = jacobian(f,u);
            E = jacobian(f,d);
            
            C = jacobian(g,x);
            D = jacobian(g,u);
            H = jacobian(g,d);
            
            f0 = subs(obj.f_sym, vertcat(x,u,d), zeros(obj.Nx + obj.Nu + obj.Nd,1));
            g0 = subs(obj.g_sym, vertcat(x,u,d), zeros(obj.Nx + obj.Nu + obj.Nd,1));
            
            %lm = LinearModel(A,B,E,C,D,H);
            lm = LinearModel();
            lm.Name = obj.Name + "_Linearized";
            copyModelProps(obj, lm);
            
            [lm.A_sym, lm.B_sym, lm.E_sym, lm.C_sym, lm.D_sym, lm.H_sym, lm.f0_sym, lm.g0_sym]...
                = deal(A,B,E,C,D,H,f0,g0);
            
            lm.init();
        end
        
        function lm = get.LinearModel(obj)
            if isempty(obj.LinearModel)
                lm = getLinearModel(obj);
                obj.LinearModel = lm;
            else
                lm = obj.LinearModel;
            end
        end
        
        function sys_h = makeSimulinkModel(obj, sys_name)
            if nargin == 1
                sys_name = 'Model_Simulink';
            end
            
            if isempty(obj.Name)
                model_name = 'Model';
                obj_name = 'Model_Object';
            else
                model_name = obj.Name;
                obj_name = sprintf('%s_Object', obj.Name);
            end
            
            try
                sys_h = load_system(sys_name);
            catch
                sys_h = new_system(sys_name);
                load_system(sys_h);
            end
            mdlWks = get_param(sys_name,'ModelWorkspace');
            assignin(mdlWks, obj_name, obj); % Assign object to Model Workspace so we can reference it
            
            % Get
            h = Simulink.findBlocks(sys_name, 'Name', model_name, 'BlockType', 'SubSystem');
            if numel(h) == 1
                model_path = getfullname(h);
            elseif numel(h) > 1
                error("Multiple subsystem blocks with name %s found in %s", model_name, sys_name);
            elseif numel(h) == 0
                % Assign model path
                model_path = sprintf('%s/%s', sys_name, model_name);
                % Copy model template into model
                add_block('Model_SimulinkTemplate/Model',model_path);
            end
            
            set_param(model_path, 'x_0', sprintf('%s.x0', obj_name));
            
            set_param([model_path,'/Model_CalcF'], 'MATLABFcn', sprintf('Model_SimulinkInterpretedFunction(u,''%s'',''%s'',@CalcFMux)', sys_name, obj_name));
            set_param([model_path,'/Model_CalcF'], 'OutputDimensions', sprintf('%s.Nx', obj_name));
            set_param([model_path,'/Model_CalcG'], 'MATLABFcn', sprintf('Model_SimulinkInterpretedFunction(u,''%s'',''%s'',@CalcGMux)', sys_name, obj_name));
            set_param([model_path,'/Model_CalcG'], 'OutputDimensions', sprintf('%s.Ny', obj_name));
            
            if obj.Nu
                set_param([model_path,'/Input1'], 'PortDimensions', sprintf('[%s.Nu,1]', obj_name));
            end
            if obj.Nd
                set_param([model_path,'/Input2'], 'PortDimensions', sprintf('[%s.Nd,1]', obj_name));
            end
        end
        
        function t = get.StateTable(obj)
            state_syms = arrayfun(@(x) string(x), obj.SymVars.x);
            t = table(state_syms, obj.StateDescriptions, 'VariableNames', ["State Variable", "Description"]);
        end
        
        function t = get.InputTable(obj)
            in_syms = arrayfun(@(x) string(x), obj.SymVars.u);
            t = table(in_syms, obj.InputDescriptions,'VariableNames', ["Input Variable", "Description"]);
        end
        
        function t = get.DisturbanceTable(obj)
            dist_syms = arrayfun(@(x) string(x), obj.SymVars.d);
            t = table(dist_syms, obj.DisturbanceDescriptions, 'VariableNames', ["Disturbance Variable", "Description"]);
        end
        
        function t = get.OutputTable(obj)
            out_syms = arrayfun(@(x) sprintf("y%d", x), 1:obj.Ny);
            t = table(out_syms', obj.OutputDescriptions, 'VariableNames', ["Output Variable", "Description"]);
        end
        
        function c = parseMuxArg(obj, mux_arg)
            c = cell(1,3);
            start_indices = 1 + cumsum([0 obj.Nx obj.Nu]);
            end_indices = cumsum([obj.Nx obj.Nu obj.Nd]);
            c{1} = mux_arg(start_indices(1):end_indices(1));
            c{2} = mux_arg(start_indices(2):end_indices(2));
            c{3} = mux_arg(start_indices(3):end_indices(3));
        end
        
        function v = CalcFMux(obj, mux_arg)
            c = parseMuxArg(obj, mux_arg);
            v = CalcF(obj, c{:});
        end
        
        function v = CalcGMux(obj, mux_arg)
            c = parseMuxArg(obj, mux_arg);
            v = CalcG(obj, c{:});
        end
        
        function exportMatlabFunctions(obj)
            f = obj.f_sym;
            g = obj.g_sym;
            vars = {[obj.SymVars.x], [obj.SymVars.u], [obj.SymVars.d]};
            if ~isempty(obj.Params)
                vars{end+1} = tunableSyms(obj.Params);
            end
            
            name = obj.Name;
            matlabFunction(f,'Vars',vars, 'File', name+"_f_sym", 'Optimize', false);
            matlabFunction(g,'Vars',vars, 'File', name+"_g_sym", 'Optimize', false);
        end
        
        function exportPythonModel(obj, opts)
            arguments
                obj
                opts.Path = pwd
                opts.Optimize = false
                opts.Name = obj.Name
            end
            Name = opts.Name;
            disp(['Exporting model ',Name,' to python...'])
            path = opts.Path;
            comp = Name;
            
            % create folder
            if not(isfolder(comp))
                mkdir(comp)
            end
            
            funcs = {obj.f_sym, obj.g_sym};
            funcs_desc = ["f", "g"];
            
            vars = {[obj.SymVars.x], [obj.SymVars.u], [obj.SymVars.d]};
            vars_desc = ["x", "u", "d"];
            param_flag = ~isempty(obj.Params);
            if param_flag
                [vars{end+1}, params] = tunableSyms(obj.Params);
                vars_desc(end+1) = "theta";
                param_tbl = dispTable(params, ["SymID", "Value", "Unit", "Description"]);
                Ntheta = numel(params);
            end
            saveLoc = path+"\"+comp+"\";
            
            % save python files
            for i = 1:numel(funcs_desc)
                fprintf('Working on %s ...\n', funcs_desc(i));
                write_py("Model_"+funcs_desc(i),"Calc_"+funcs_desc(i),funcs{i}, vars, vars_desc)
            end
            
            % save Jacobians
            J = cell(numel(funcs_desc), numel(vars_desc));
            JacStruct = struct();
            for i = 1:numel(funcs_desc)
                for j = 1:numel(vars_desc)
                    jac_name = funcs_desc(i)+"_"+vars_desc(j);
                    fprintf("Writing Jacobian %s\n", jac_name);
                    
                    J{i,j} = jacobian(funcs{i}, vars{j});
                    [J_flat, st] = flattenJac(jac_name, J{i,j});
                    JacStruct.(funcs_desc(i)).(vars_desc(j)) = st; % Create list of structures with information about Jacobian structure.  
                    if st.NCalc > 0
                        write_py("ModelJ_"+jac_name,"CalcJ_"+jac_name,J_flat, vars, vars_desc);
                    end
                end
            end
            
            % Save Model Information
            metadata_struct = struct();
            metadata_struct.Name = obj.Name;
            
            metadata_struct.Nx = obj.Nx;
            metadata_struct.Nu = obj.Nu;
            metadata_struct.Nd = obj.Nd;
            metadata_struct.Ny = obj.Ny;
            if param_flag
                metadata_struct.Ntheta = Ntheta;
            end
            
            metadata_struct.x0 = obj.x0;
            
            metadata_struct.StateTable = table2struct(obj.StateTable);
            metadata_struct.InputTable = table2struct(obj.InputTable);
            metadata_struct.DisturbanceTable = table2struct(obj.DisturbanceTable);
            metadata_struct.OutputTable = table2struct(obj.OutputTable);
            if param_flag
                param_tbl = dispTable(params, ["SymID", "Value", "Unit", "Description"]);
                metadata_struct.ParamTable = table2struct(param_tbl);
            end
            metadata_struct.JacStruct = JacStruct;
            
            metadata_json = jsonencode(metadata_struct,'PrettyPrint',true);
            fid=fopen(saveLoc+"ModelMetadata.json",'w');
            fprintf(fid, metadata_json);
            fclose(fid);
            
            function [J_flat, st] = flattenJac(Name,J)
                % this function accounts for a significant amount of the function call
                % time. Maybe figure out how to simplify this
                idx = hasSymType(J,'variable');
                [rowsCalc, colsCalc] = find(idx);
                idxCalc = find(idx);
                [rowsConst, colsConst] = find(~idx);
                [idxCon] = find(~idx);
                valCon = double(J(idxCon)); % this line accounts for like 90% of the computation time of this function. Would be nice to find an alternative
                idxZero = valCon==0;
                rowsConst(idxZero) = [];
                colsConst(idxZero) = [];
                valConst = valCon(~idxZero);
                
                st = struct();
                st.Name = Name;
                st.NCalc = numel(idxCalc);
                st.NConst = numel(valConst);
                st.iCalc = idxCalc;
                st.rCalc = rowsCalc;
                st.cCalc = colsCalc;
                st.rConst = rowsConst;
                st.cConst = colsConst;
                st.valConst = valConst;
                
                J_flat = J(idx);
            end
            
            function write_py(filename,functname,symFunc, VarCell, VarDesc) % May also need nn here
                filepath = string(saveLoc)+"\"+filename;
                
                % this block turns the symbolic expression into text that can be used to
                % generate the ML function
                n = length(symFunc);
                out = sym('out',[n,1]);
                func = out == symFunc;
                strFunc = sprintf('\t%s \n',func);
                oName = "out"+string(1:n);
                
                %%
                % Replace Header
                % CalcText(1:8) = [];
                funcDef = ["import math"];
                funcDef(end+1,1) = ["import numpy as np"];
                funcDef(end+1,1) = [""];
                funcDef(end+1,1) = ["def " + functname + "(" + strjoin(VarDesc, ",") + "):"]; % Function arguments are x,u,d,theta
                funcDef(end+1,1) = "# auto-generated function from matlab";
                funcDef(end+1,1) = [""];
                
                % replace input variable parsing
                func_vars = symvar(symFunc);
                ParseText = string.empty();
                for m = 1:numel(VarCell)
                    sym_vars = VarCell{m}; 
                    str_vars = string(sym_vars);
                    for n = 1:numel(sym_vars)
                        if ismember(sym_vars(n), func_vars)
                            index_str = VarDesc(m) + "[" + string(n - 1) + "]";
                            eq_str = sprintf("%s = %s", str_vars(n), index_str); 
                            ParseText = [ParseText, eq_str];
                        end
                    end
                end
                
                % replace matlab syntax
                strFunc = strrep(strFunc,'==','=');
                strFunc = strrep(strFunc,'*','*');
                strFunc = strrep(strFunc,'/','/');
                strFunc = strrep(strFunc,'^','**');
                strFunc = strrep(strFunc,'sqrt','math.sqrt');
                strFunc = strrep(strFunc,'pi','np.pi');
                
                % function return string
                returnFunc = ['return ', sprintf('%s, ',oName{1:end-1}),oName{end}];
                
                % write to file
                fid = fopen(filepath + ".py",'w');
                fprintf(fid,'%s\n',funcDef);
                fprintf(fid,'\t%s\n',ParseText);
                fprintf(fid,'%s',strFunc);
                fprintf(fid,'\t%s\n',returnFunc);
                fclose(fid);
            end
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
            
            cell_flag = isa(syms, 'cell');
            if cell_flag
                funcs = cellfun(@processSym, syms, 'UniformOutput', false);
            else
                funcs = processSym(syms);
            end
            
            function func = processSym(sym)
                if isa(sym, 'sym')
                    if ~isempty(obj.Params)
                        func = matlabFunction(obj.Params, sym, vars);
                    else
                        func = matlabFunction(sym,'Vars',vars);
                    end
                else
                    func = @(varargin) sym;
                end
            end
        end
        
        function copyModelProps(obj_from, obj_to, opts)
            arguments
                obj_from
                obj_to
                opts.Properties = ["Nx","Nu","Nd","Ny","x0",...
                    "StateDescriptions", "InputDescriptions", "DisturbanceDescriptions", "OutputDescriptions",...
                    "SymVars","Params"];
            end
            
            for prop = opts.Properties
                val = obj_from.(prop);
                if ~isempty(val)
                    obj_to.(prop) = val;
                end
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

