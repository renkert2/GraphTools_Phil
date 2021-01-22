classdef GraphModel < Model
    %GRAPHMODEL Contains all of the operations useful for working with the
    %graph models.  Most of the code will go here.  
    %   Detailed explanation goes here
    properties
        graph Graph = Graph.empty()
    end
    
    properties (SetAccess = private)
        C_coeff % capacitance coefficient matrix
        CType Type_Capacitance = Type_Capacitance.empty()
        P_coeff % capacitance coefficient matrix
        PType Type_PowerFlow = Type_PowerFlow.empty()
        CapFunction (:,1) LookupFunction = LookupFunction.empty()
        x_init % capacitance coefficient matrix
        DynType DynamicTypes = DynamicTypes.EnergyFlow
        D % capacitance coefficient matrix
        B (:,:,:)% input mapping matrix   
    end
    
    methods
        function obj = GraphModel(varargin)
            if nargin == 1
                obj.graph = varargin{1};
                init(obj);
            end
        end      

        function init(obj)
            % make vertex matrices
            obj.x_init  = vertcat(obj.graph.Vertices.Initial); 
            obj.DynType = vertcat(obj.graph.Vertices.DynamicType); 
            
            % D matrix
            Dmat = zeros(obj.graph.v_tot,obj.graph.Nee);
            E_idx = arrayfun(@(x) find(x==obj.graph.InternalVertices),vertcat(obj.graph.ExternalEdges.HeadVertex));            
            for i  = 1:length(E_idx)
                Dmat(E_idx(i),i) = 1;
            end
            obj.D = Dmat;
            
            % B matrix
            Eint = obj.graph.InternalEdges;
            numU = max(arrayfun(@(x) length(x.Input),Eint));

            obj.B = zeros(obj.graph.Ne, obj.graph.Nu, numU); % Inputs added along second dimensions; separate inputs along third dimension
            for i = 1:numel(Eint)
                for j = 1:numel(Eint(i).Input)
                    obj.B(i,:,j) = (Eint(i).Input(j)==obj.graph.Inputs);
                end
            end
            
            % C matrix
            CTypeAll = vertcat(obj.graph.Vertices(:).Capacitance);
            numCType = arrayfun(@(x) length(x.Capacitance),obj.graph.Vertices);
            [obj.C_coeff,obj.CType] = MakeCoeffMatrix(obj.graph.Vertices,CTypeAll,numCType);

            % P matrix
            PTypeAll = vertcat(Eint(:).PowerFlow); % list of all capacitance types
            numPType = arrayfun(@(x) length(x.PowerFlow),Eint); % find number of capacitance types per vertex
            [obj.P_coeff,obj.PType] = MakeCoeffMatrix(Eint,PTypeAll,numPType);
            
            % Lookup Functions
            obj.CapFunction = vertcat(obj.graph.InternalVertices(:).CapFunction);
%             Functions = vertcat(obj.graph.InternalVertices(:).CapFunction);
%             obj.CapFunction = Functions(~arrayfun(@(x) isempty(x.Breakpoints),Functions));
         
            obj.Nx = sum(any(obj.C_coeff ~= 0,2));
            obj.Nu = obj.graph.Nu;
            obj.Nd = obj.graph.Nev + obj.graph.Nee;
            obj.SymbolicSolve
            
            init@Model(obj);
            
            
        end
        
        function Modify(obj)
        
        end
        
        function Simulate(obj)
            
        end
        
        function x = StateNames(obj)
            Desc = vertcat(obj.graph.DynamicVertices.Description);
            Blks = vertcat(vertcat(obj.graph.DynamicVertices.Parent).Name);
            x = join([Blks,repmat('\',length(Blks),1),Desc]);
        end
        
        function x = InputNames(obj)
            Desc = vertcat(obj.graph.Inputs.Description);
            Blks = vertcat(vertcat(obj.graph.Inputs.Parent).Name);
            x = join([Blks,repmat('\',length(Blks),1),Desc]);
        end
        
        function x = DisturbanceNames(obj)
            Desc = [vertcat(obj.graph.ExternalVertices.Description);vertcat(obj.graph.ExternalEdges.Description)];
            Blks = [vertcat(vertcat(obj.graph.ExternalVertices.Parent).Name); vertcat(vertcat(obj.graph.ExternalEdges.Parent).Name)];
            x = join([Blks,repmat('\',length(Blks),1),Desc]);
        end
        
        function x = OutputNames(obj)
            Desc = vertcat(obj.graph.InternalVertices.Description);
            Blks = vertcat(vertcat(obj.graph.InternalVertices.Parent).Name);
            x = join([Blks,repmat('\',length(Blks),1),Desc]);
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
            
            h = plot(obj.graph,varargin{:});
            
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
                labelnode(h,1:obj.graph.Nv,obj.OutputNames)
            end

            function LabelEdges(obj,h)
                % Get state vector filled up with symbolic varialbes
                x       = sym('x%d'      ,[obj.graph.v_tot        1]); % dynamic states
                u       = sym('u%d'      ,[obj.graph.Nu           1]); % inputs
                % Calculate power flows and capacitances
                P = CalcP(obj,x,u); % calculates power flows
                labeledge(h,obj.graph.E(:,1)',obj.graph.E(:,2)',string(P)) %label Edges
            end
            
            function LabelDisturbances(obj,h)
                labelnode(h,obj.graph.Nv+1:obj.graph.v_tot,obj.DisturbanceNames(1:obj.graph.Nev))
                labeledge(h,obj.graph.Ne+1:obj.graph.Ne+obj.graph.Nee,obj.DisturbanceNames(obj.graph.Nev+1:end))  
            end
        
        end

          
        function SymbolicSolve(obj) % this function will only work for symbolic expressions at the moment
        
            idx_x_d = (sum(abs(obj.C_coeff(1:obj.graph.Nv,:)),2) ~= 0);
            idx_x_a = (sum(abs(obj.C_coeff(1:obj.graph.Nv,:)),2) == 0);
            idx_x_e = obj.graph.Nv+1:obj.graph.Nv+obj.graph.Nev;
            
            x       = sym('x%d'      ,[sum(idx_x_d)        1]); % dynamic states
            x_a     = sym('x_a%d'    ,[sum(idx_x_a)        1]); % algebraic states
            u       = sym('u%d'      ,[obj.graph.Nu              1]); % inputs
            d       = sym('d%d'      ,[obj.graph.Nev+obj.graph.Nee 1]);
            x_e     = d(1:obj.graph.Nev); % external states
            P_e     = d(obj.graph.Nev+1:end);
            
            % Get state vector filled up with symbolic varialbes
            x_full(idx_x_d,1) = x;
            x_full(idx_x_a,1) = x_a;
            x_full(idx_x_e,1) = x_e;
            
            % Calculate power flows and capacitances
            P = CalcP(obj,x_full,u); % calculates power flows
            C = CalcC(obj,x_full); % calcualtes capacitance
            
            % Solve system dynamics           
            eqnA(1:sum(idx_x_a),1) = -obj.graph.M(idx_x_a,:)*P + obj.D(idx_x_a,:)*P_e == 0; % system of algebraic equations
            eqnD(1:sum(idx_x_d),1) = diag(C(idx_x_d))^-1*(-obj.graph.M(idx_x_d,:)*P + obj.D(idx_x_d,:)*P_e); % system of dynamic equations (
            [A,Bu] = equationsToMatrix(eqnA,x_a); % convert eqnA to the form Ax=B
            x_a_solution = linsolve(A,Bu); % find solution to the algebraic system
            x_d_solution = subs(eqnD,x_a,x_a_solution); % plug in the algebraic system solution into the dynamic system equations
                        
            % Store symbolic calculations
            obj.f_sym = x_d_solution; % system derivatives
            obj.g_sym = [x_full(idx_x_d);x_a_solution]; % all system states

            
        end
            
        function [P] = CalcP(obj,x0,u0)
            % CalcP calculates the power flows of a graph model.
            
            %%% INPUTS
            % Sys  - System graph model object
            % x0   - state vector
            % u0   - input vector
            % Pmap - vector of lookup map values
            
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
            
            % P   = zeros(size(Sys.P_coeff_mod));
            xt = obj.graph.Tails*x0; %tail states
            xh = obj.graph.Heads*x0; %head states
            [~,~,Nu] = size(obj.B); % max number of inputs incident per edge
            Ne = obj.graph.Ne;
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
            
        end
        
        function [C] = CalcC(obj,x0)
            % CalcC calculates the capacitance values of a graph model.
            
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
                c(:,i) = obj.C_coeff(:,i).*obj.CType(i).calcVal(x0); % the 1 and 0 in these lines will need to be changed
            end
            c = sum(c,2); % sum across capacitance coefficients
            
            % function lookups
            LF = sym(ones(obj.graph.Nv,1));
            for i = 1:length(LF)
                if ~isempty(obj.CapFunction(i).Breakpoints)
                    [~,idx] = ismember(obj.CapFunction(i).Breakpoints(:),obj.graph.Vertices);
                    input = num2cell(x0(idx));
                    LF(i) = obj.CapFunction(i).Function.calcVal(input{:});
                end
            end
            
            
            C = LF.*c(1:obj.graph.Nv); % solve for the capacitance of each vertex
            
        end
        
    end
        
end

