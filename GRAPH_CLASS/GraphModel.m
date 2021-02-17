classdef GraphModel < Model
    %GRAPHMODEL Contains all of the operations useful for working with the
    %graph models.  Most of the code will go here.  
    %   Detailed explanation goes here
    properties
        Graph Graph = Graph.empty()
        AutomaticModify logical = true
        CalcPMethod CalcPMethods = CalcPMethods.Default
    end
    
    properties (SetAccess = private)
        C_coeff % capacitance coefficient matrix
        CType Type_Capacitance = Type_Capacitance.empty()
        P_coeff % capacitance coefficient matrix
        PType Type_PowerFlow = Type_PowerFlow.empty()
        x_init % capacitance coefficient matrix
        DynType DynamicTypes = DynamicTypes.EnergyFlow
        D % capacitance coefficient matrix
        B (:,:,:)% input mapping matrix
        
        P_sym (:,1) sym % Symbolic Representation of Power Flows
    end
    
    properties (SetAccess = protected, GetAccess = protected)
        CalcP_func % Matlab function of Power Flows
    end
    
    properties (Dependent)
        VertexTable table
        EdgeTable table
    end

    methods
        function obj = GraphModel(arg1, opts)
            arguments 
                arg1 (1,1)
                opts.Linearize logical = true
                opts.CalcPMethod CalcPMethods = "Default"
                opts.SymParams_HandleMethod SymParams_HandleMethods = SymParams_HandleMethods.AugmentMatlabFunctions
            end
            
            if isa(arg1,'Graph')
                obj.Graph = arg1;
            elseif isa(arg1,'Component')
                obj.Graph = arg1.Graph;
            else
                error('Invalid argument to GraphModel.  Must be of type Graph or Component')
            end
            
            obj.CalcPMethod = opts.CalcPMethod;
            obj.LinearizeFlag = opts.Linearize;
            obj.SymParams_HandleMethod = opts.SymParams_HandleMethod;
            
            init(obj);
        end      

        function init(obj)
            init@Model(obj); 
        end
        
        function initSymbolic(obj)
            obj.Nx = obj.Graph.Nx;
            obj.Nu = obj.Graph.Nu;
            obj.Nd = obj.Graph.Nev + obj.Graph.Nee;
            obj.Ny = obj.Graph.Nv + numel(obj.Graph.Outputs);
            
            obj.SymParams = obj.Graph.SymParams;
            obj.SymParams_Vals = obj.Graph.SymParams_Vals;
            
            % make vertex matrices
            if ~isempty(obj.Graph.DynamicVertices)
                obj.x_init = vertcat(obj.Graph.DynamicVertices.Initial);
            else
                obj.x_init = [];
            end
 
            obj.DynType = vertcat(obj.Graph.Vertices.DynamicType); 
            
            % D matrix
            Dmat = zeros(obj.Graph.v_tot,obj.Graph.Nee);
            E_idx = arrayfun(@(x) find(x==obj.Graph.InternalVertices),vertcat(obj.Graph.ExternalEdges.HeadVertex));            
            for i  = 1:length(E_idx)
                Dmat(E_idx(i),i) = 1;
            end
            obj.D = Dmat;
            
            % B matrix
            Eint = obj.Graph.InternalEdges;
            numU = max(arrayfun(@(x) length(x.Input),Eint));

            obj.B = zeros(obj.Graph.Ne, obj.Graph.Nu, numU); % Inputs added along second dimensions; separate inputs along third dimension
            for i = 1:numel(Eint)
                for j = 1:numel(Eint(i).Input)
                    obj.B(i,:,j) = (Eint(i).Input(j)==obj.Graph.Inputs);
                end
            end
            
            % C matrix
            % Merge Note - May need to change from x.InternalVertices to x.Vertices
            CTypeAll = vertcat(obj.Graph.InternalVertices(:).Capacitance);
            numCType = arrayfun(@(x) length(x.Capacitance),obj.Graph.InternalVertices);
            [obj.C_coeff,obj.CType] = MakeCoeffMatrix(obj.Graph.InternalVertices,CTypeAll,numCType);

            % P matrix
            PTypeAll = vertcat(Eint(:).PowerFlow); % list of all capacitance types
            numPType = arrayfun(@(x) length(x.PowerFlow),Eint); % find number of capacitance types per vertex
            [obj.P_coeff,obj.PType] = MakeCoeffMatrix(Eint,PTypeAll,numPType);
            
            initSymbolic@Model(obj);
            
            obj.SymbolicSolve;
        end
        
        function initNumerical(obj)                
            u_mod = genSymVars('u%d',max([obj.Graph.Nu,2])); % inputs - modified from SymVars.u to force MATLAB to vectorize CalcP_func even if there's a single input
            vars = {[obj.SymVars.x_full], [u_mod]};
            obj.CalcP_func = genMatlabFunctions(obj, obj.P_sym, vars);
            
            initNumerical@Model(obj);
        end
        
        function x = defineStateNames(obj)
            if ~isempty(obj.Graph.DynamicVertices)
                Desc = vertcat(obj.Graph.DynamicVertices.Description);
                Blks = vertcat(vertcat(obj.Graph.DynamicVertices.Parent).Name);
                x = join([Blks,repmat('\',length(Blks),1),Desc]);
            else
                x = string.empty();
            end
            
        end
        
        function x = defineInputNames(obj)
            if ~isempty(obj.Graph.Inputs)
                Desc = vertcat(obj.Graph.Inputs.Description);
                Blks = vertcat(vertcat(obj.Graph.Inputs.Parent).Name);
                x = join([Blks,repmat('\',length(Blks),1),Desc]);
            else
                x = string.empty();
            end
        end
        
        function x = defineDisturbanceNames(obj)
            if ~isempty(obj.Graph.ExternalVertices)
                ext_verts_desc = vertcat(obj.Graph.ExternalVertices.Description);
                ext_verts_parents = vertcat(vertcat(obj.Graph.ExternalVertices.Parent).Name);
            else
                ext_verts_desc = [];
                ext_verts_parents = [];
            end
            
            if ~isempty(obj.Graph.ExternalEdges)
                ext_edges_desc = vertcat(obj.Graph.ExternalEdges.Description);
                ext_edges_parents =  vertcat(vertcat(obj.Graph.ExternalEdges.Parent).Name);
            else
                ext_edges_desc = [];
                ext_edges_parents = [];
            end
                
            Desc = [ext_verts_desc; ext_edges_desc];
            Blks = [ext_verts_parents; ext_edges_parents];
            x = join([Blks,repmat('\',length(Blks),1),Desc]);
        end
        
        function x = defineOutputNames(obj)
            DescX = vertcat(obj.Graph.InternalVertices.Description);
            BlksX = vertcat(vertcat(obj.Graph.InternalVertices.Parent).Name);
            
            if ~isempty(obj.Graph.Outputs)
                DescY = vertcat(obj.Graph.Outputs.Description);
                BlksY = vertcat(vertcat(obj.Graph.Outputs.Parent).Name);
            end
            Desc = [DescX; DescY];
            Blks = [BlksX; BlksY];
            x = join([Blks,repmat('\',length(Blks),1),Desc]);
        end
        
        function vertex_table = get.VertexTable(obj)
            names = vertcat(obj.Graph.Vertices.Description);
            parents = vertcat(vertcat(obj.Graph.Vertices.Parent).Name);
            types = vertcat(obj.Graph.Vertices.VertexType);
            vertex_table = table((1:(obj.Graph.Nv+obj.Graph.Nev))', parents, names, types, 'VariableNames', ["Vertices", "Component", "Description", "Domain"]);
        end
        
        function edge_table = get.EdgeTable(obj)
            digits(4)
            pflows = vpa(obj.P_sym);
            pflows_strings = arrayfun(@string, pflows);
            digits(32)
            
            parents = vertcat(vertcat(obj.Graph.InternalEdges.Parent).Name);
            
            edge_table = table((1:(obj.Graph.Ne - obj.Graph.Nee))',parents, pflows_strings, 'VariableNames', ["Edges", "Component", "PowerFlows"]);
        end
        
        function [t,x, pf] = Simulate(obj, inputs, disturbances, params, t_range, opts)
            % SIMULATE(GraphModel, inputs, disturbances, t_range, opts)
            % inputs and disturbances must be column vectors of appropriate size.
            % Inputs and disturbances can be anonymous functions of time or constant values
            % If the Graph includes symbolic parameters, pass numerical values to the params argument as a column vector.  Leave as an empty double [] if no symbolic parameters exist
            % Simulate usees the first and last entries of t_range if dynamic states are calculated
            % with ODC23t, or the entire t_range vector if only algebraic states are calculated
            
            arguments
                obj
                inputs
                disturbances
                params
                t_range
                opts.PlotStates logical = true
                opts.PlotInputs logical = false
                opts.PlotDisturbances logical = false
                opts.StateSelect = []
                opts.Solver = @ode23t
                opts.SolverOpts struct = struct.empty()
            end
            
            input_function_flag = isa(inputs, 'function_handle');
            disturbance_function_flag = isa(disturbances, 'function_handle');
            
            if ~isempty(obj.Graph.DynamicVertices)
                % Dynamic states exist
                xdot = processArgs(@CalcF, inputs, disturbances, params); % Might need to change processArgs so things run faster
                [t,xdyn] = opts.Solver(xdot, [t_range(1) t_range(end)], obj.x_init, opts.SolverOpts);
                if ~isempty(obj.Graph.AlgebraicVertices)
                    xfull = processArgs(@CalcG, inputs, disturbances, params);
                    x = xfull(t',xdyn')';
                else
                    x = xdyn;
                end
            else
                % Internal States are all constant
                xfull = processArgs(@CalcG, inputs, disturbances);
                t = t_range;
                x = xfull(t,[])';    
            end
            
            if nargout == 3
                if input_function_flag
                    pf = CalcP(obj, x', inputs(t)', repmat(params,1,numel(t)))';
                else
                    pf = CalcP(obj, x', repmat(inputs,1,numel(t)), repmat(params,1,numel(t)))';
                end
            end
            
            if any([opts.PlotStates opts.PlotInputs opts.PlotDisturbances])
                hold on
                lgnd = string.empty();
                if opts.PlotStates
                    if opts.StateSelect
                        x = x(:,opts.StateSelect);
                        names = obj.StateNames(opts.StateSelect);
                    else
                        names = obj.StateNames;
                    end
                    plot(t,x)
                    lgnd = vertcat(lgnd, names);
                end
                
                if opts.PlotInputs && ~isempty(inputs)
                    if input_function_flag
                        plot(t,inputs(t))
                    else
                        plot(t,inputs*ones(size(t)));
                    end
                    lgnd = vertcat(lgnd,obj.InputNames);
                end   
                
                if opts.PlotDisturbances && ~isempty(disturbances)
                    if disturbance_function_flag
                        plot(t,disturbances(t')')
                    else
                        plot(t,(disturbances*ones(size(t))')');
                    end
                    lgnd = vertcat(lgnd,obj.DisturbanceNames);
                end
                legend(lgnd)
                hold off
            end
            
            function xfunc = processArgs(func, inputs_arg, disturbances_arg, params_arg)               
                if input_function_flag && disturbance_function_flag
                    xfunc = @(t,x) func(obj, x, inputs_arg(t), disturbances_arg(t),repmat(params_arg,1,numel(t)));
                elseif input_function_flag
                    xfunc = @(t,x) func(obj, x, inputs_arg(t), repmat(disturbances_arg,1,numel(t)),repmat(params_arg,1,numel(t)));
                elseif disturbance_function_flag
                    xfunc = @(t,x) func(obj, x, repmat(inputs_arg,1,numel(t)), disturbances_arg(t),repmat(params_arg,1,numel(t)));
                else
                    xfunc = @(t,x) func(obj, x, repmat(inputs_arg,1,numel(t)), repmat(disturbances_arg,1,numel(t)),repmat(params_arg,1,numel(t)));
                end
            end
        end

        function h = plot(obj,varargin)
            % Pass the Graph Model you would like to plot with modifiers in
            % NAME-VALUE pairs. The modifiers can include ANY modifiers
            % used in the MATLAB digraph plotting PLUS the added option for
            % graph model detailed labels. To invoke the detailed labels 
            % modifiers, use the NAME "DetailedLabels" with VALUE,
            % "States","Edges","Disturbances", or "All".
            
            if mod(length(varargin),2) == 1
                error('Plot modifiers must be in name-value pairs')
            end      
            % extract and remove the DetailedLabels options from varargin 
            idxName = find(cellfun(@(x) strcmp('DetailedLabels',x),varargin));
            idxValue = idxName+1;
            LABELS = {varargin{idxValue}};
            varargin([idxName,idxValue]) = [];
            
            h = plot(obj.Graph,varargin{:});
            
            for i = 1:length(LABELS)
                opts = {'All','States','Edges','Disturbances'};             
                switch LABELS{i}
                    case opts{1}
                        LabelStates(obj,h)                    
                        LabelEdges(obj,h)                    
                        LabelDisturbances(obj,h)                    
                    case opts{2}
                        LabelStates(obj,h)                    
                    case opts{3}
                        LabelEdges(obj,h)                    
                    case opts{4}                 
                        LabelDisturbances(obj,h)
                    otherwise
                        error([sprintf('Provide a valid argument value for DetailedLabels:\n') sprintf('-%s\n',opts{:})]) % update this to list valid arguments
                end
            end
            
            set(gcf,'WindowButtonDownFcn',@(f,~)edit_graph(f,h))
            
            function LabelStates(obj,h)
                labelnode(h,1:obj.Graph.Nv,obj.OutputNames(1:obj.Graph.Nv))
            end

            function LabelEdges(obj,h)
                % Get state vector filled up with symbolic varialbes
                x       = sym('x%d'      ,[obj.Graph.v_tot        1]); % dynamic states
                u       = sym('u%d'      ,[obj.Graph.Nu           1]); % inputs
                % Calculate power flows and capacitances
                P = CalcP_Sym(obj,x,u); % calculates power flows
                labeledge(h,obj.Graph.E(:,1)',obj.Graph.E(:,2)',string(P)) %label Edges
            end
            
            function LabelDisturbances(obj,h)
                labelnode(h,obj.Graph.Nv+1:obj.Graph.v_tot,obj.DisturbanceNames(1:obj.Graph.Nev))
                labeledge(h,obj.Graph.Ne+1:obj.Graph.Ne+obj.Graph.Nee,obj.DisturbanceNames(obj.Graph.Nev+1:end))  
            end
        end
          
        function SymbolicSolve(obj) % this function will only work for symbolic expressions at the moment
            if ~isa(obj.C_coeff, 'sym')
                idx_x_d = (sum(abs(obj.C_coeff(1:obj.Graph.Nv,:)),2) ~= 0);
                idx_x_a = ~idx_x_d;
            else
                C_sum = sum(abs(obj.C_coeff(1:obj.Graph.Nv,:)),2);
                idx_x_a = arrayfun(@(x) isequal(x,sym(0)), C_sum);
                idx_x_d = ~idx_x_a;
            end
            %idx_x_e = obj.Graph.Nv+1:obj.Graph.Nv+obj.Graph.Nev;
 
            x = obj.SymVars.x; % Dynamic States 
            x_a = genSymVars('x_a%d', sum(idx_x_a)); % Algebraic States
            u = obj.SymVars.u; % Inputs
            d = obj.SymVars.d; % Disturbances
            x_e     = d(1:obj.Graph.Nev); % external states
            
            P_e     = d(obj.Graph.Nev+1:end);
            
            %Get state vector filled up with symbolic varialbes
            if any(idx_x_d)
                x_internal(idx_x_d,1) = x;
            end
            if any(idx_x_a)
                x_internal(idx_x_a,1) = x_a;
            end
            
            x_full = vertcat(x_internal, x_e);
            
            obj.SymVars.x_full = x_full; % Add full list of symbolic state variables, used later in initNumerical()
            
            % Calculate power flows and capacitances
            P = CalcP_Sym(obj,x_full,u); % calculates power flows
            obj.P_sym = P;
            
            C = CalcC_Sym(obj,x_full, x_internal); % calcualtes capacitance
            
            % Solve system dynamics
            % Process Algebraic States
            if any(idx_x_a)
                eqnA_temp = -obj.Graph.M(idx_x_a,:)*P + obj.D(idx_x_a,:)*P_e;
                if obj.AutomaticModify
                    for i = 1:length(eqnA_temp)
                        eqn = eqnA_temp(i);
                        factors = factor(eqn);
                        if ismember(x_a(i), factors)
                            eqn = simplify(eqn/x_a(i));
                        end
                        eqnA_temp(i,1) = eqn;
                    end
                end
                eqnA(1:sum(idx_x_a),1) = eqnA_temp == 0; % system of algebraic equations
                [A,Bu] = equationsToMatrix(eqnA,x_a); % convert eqnA to the form Ax=B
                x_a_solution = linsolve(A,Bu); % find solution to the algebraic system
            else
                x_a_solution = sym.empty();
            end
            
            if any(idx_x_d)
                eqnD(1:sum(idx_x_d),1) = diag(C(idx_x_d))^-1*(-obj.Graph.M(idx_x_d,:)*P + obj.D(idx_x_d,:)*P_e); % system of dynamic equations (
                if any(idx_x_a)
                    x_d_solution = simplifyFraction(subs(eqnD,x_a,x_a_solution)); % plug in the algebraic system solution into the dynamic system equations
                else
                    x_d_solution = simplifyFraction(eqnD);
                end
            else
                x_d_solution = sym.empty();
            end
            
            if ~isempty(obj.Graph.Outputs)
                Y = CalcY_Sym(obj,x_full,u);
                if any(idx_x_a)
                    Y = simplifyFraction(subs(Y,x_a,x_a_solution)); % plug in the algebraic system solution
                else
                    Y = simplifyFraction(Y);
                end   
            else
                Y = [];
            end
                        
            % Store symbolic calculations
            obj.f_sym = x_d_solution; % system derivatives
            obj.g_sym = [x_full(idx_x_d);x_a_solution;Y]; % all system states       
        end
            
        function [P] = CalcP_Sym(obj,x0,u0)
            % CalcP_Sym calculates the power flows of a Graph model.
            
            %%% INPUTS
            % Sys  - System Graph model object
            % x0   - state vector
            % u0   - input vector
            % opts.
            % - opts.Method = Default: vectorized calculation of every power flow
            % - opts.Method = Edges: calculates powerflows independently for each edge.  Avoids some numerical issues with complicated powerflows.  
            
            %%% OUTPUTS
            % P - Power flows
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Author: Christopher T. Aksland
            % Association: University of Illionis at Urbana-Champaign
            % Contact: aksland2@illinois.edu
            % Revision History:
            % 11/22/2020 - Function creation
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Potential improvements
            % -
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            if obj.CalcPMethod == "Default"
                xt = obj.Graph.Tails*x0; %tail states
                xh = obj.Graph.Heads*x0; %head states
                [~,~,Nu] = size(obj.B); % max number of inputs incident per edge
                Ne = obj.Graph.Ne;
                u = sym(zeros(Ne,Nu)); % initialize edge input data.
                for i = 1:Nu
                    u(:,i) = obj.B(:,:,i)*u0;
                end

                % calculate the powerflow along each edge. Note the 3x vector size from
                % repmat required to simulate a multi-domain graph
                P = sym(zeros(size(obj.P_coeff)));
                for i = 1:size(obj.P_coeff,2)
                    P(:,i) = obj.P_coeff(:,i).*obj.PType(i).calcVal(xt,xh,u);
                end

                % sum the powerflow coefficients
                P = sum(P,2);
            elseif obj.CalcPMethod == "Edges"
                xt = obj.Graph.Tails*x0; %tail states
                xh = obj.Graph.Heads*x0; %head states
                Ne = obj.Graph.Ne - obj.Graph.Nee;
                P = sym(zeros(Ne,1));
                for i = 1:Ne
                    edge = obj.Graph.InternalEdges(i);
                    types = edge.PowerFlow;
                    coeffs = edge.Coefficient;
                    u = squeeze(obj.B(i,:,:)).'*u0; % Column vector of inputs corresponding to this edge
                    types_sym = arrayfun(@(x) x.calcVal(xt(i),xh(i),u.'),types);
                    pflows = coeffs.'*types_sym;
                    P(i,1) = pflows;
                end      
            end
        end
        
        function [C] = CalcC_Sym(obj,x_full, x_internal)
            % CalcC_Sym calculates the capacitance values of a graph model.
            
            %%% INPUTS
            % Sys  - System graph model object
            % x0   - state vector
            
            %%% OUTPUTS
            % C - Capacitance vector
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Author: Christopher T. Aksland
            % Association: University of Illionis at Urbana-Champaign
            % Contact: aksland2@illinois.edu
            % Revision History:
            % 11/22/2020 - Function creation
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Potential improvements
            % -
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % graph lookups
            c   = sym(zeros(size(obj.C_coeff)));
            % caculate the capacitance of each vertex for each coefficient
            for i = 1:size(obj.C_coeff,2)
                c(:,i) = obj.C_coeff(:,i).*obj.CType(i).calcVal(x_internal); % the 1 and 0 in these lines will need to be changed
            end
            c = sum(c,2); % sum across capacitance coefficients
            
            % function lookups
            LF = sym(ones(obj.Graph.Nv,1));
            for i = 1:length(LF)
                cf = obj.Graph.Vertices(i).CapFunction;
                if ~isempty(cf)
                    [~,idx] = ismember(cf.Breakpoints(:),obj.Graph.Vertices);
                    input = num2cell(x_full(idx));
                    LF(i) = cf.Function.calcVal(input{:});
                end
            end
            
            
            C = LF.*c(1:obj.Graph.Nv); % solve for the capacitance of each vertex  
        end
        
        function [Y] = CalcY_Sym(obj,x0,u0)
            % CalcY_Sym calculates the outputs of a graph model.              
                    
            % function lookups
            OF = sym(ones(length(obj.Graph.Outputs),1));
            for i = 1:length(OF)
                of = obj.Graph.Outputs(i);
                input = [];
                if ~isempty(of)
                    for j = 1:length(of.Breakpoints)
                        if isa(of.Breakpoints{j},'GraphVertex')                           
                            [~,idx] = ismember(of.Breakpoints{j},obj.Graph.Vertices); 
                            input = [input x0(idx)];
                        elseif isa(of.Breakpoints{j},'GraphInput')     
                            [~,idx] = ismember(of.Breakpoints{j},obj.Graph.Inputs);
                            input = [input u0(idx)];
                        else
                            error('Breakpoint not Vertex or Input object.')
                        end
                    
                    end
                    input = num2cell(input);
                    OF(i) = of.Function.calcVal(input{:});
                end
            end
            
            
            Y = OF; % solve for the capacitance of each vertex  
        end
        
        function P = CalcP(obj, x_full, u, params)
            param_lengths = [obj.Graph.Nv, obj.Nu, numel(obj.SymParams)];
            
            if nargin == 3 || isempty(obj.SymParams)
                vars = {x_full, u};
            elseif nargin == 4
                vars = {x_full, u, params};
            end
            
            for i = 1:numel(vars)
                assert(size(vars{i},1) >= param_lengths(i), "Argument %d requires %d entries", i, param_lengths(i));
            end
            
            P = CalcX(obj, obj.CalcP_func, vars);
        end
        
    end
        
        
end

