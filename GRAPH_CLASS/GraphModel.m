classdef GraphModel < Model
    %GRAPHMODEL Contains all of the operations useful for working with the
    %graph models.  Most of the code will go here.  
    %   Detailed explanation goes here
    properties
        graph Graph = Graph.empty()
        
        C_coeff % capacitance coefficient matrix
        CType Type_Capacitance = Type_Capacitance.empty()
        P_coeff % capacitance coefficient matrix
        PType Type_PowerFlow = Type_PowerFlow.empty()
        x_init % capacitance coefficient matrix
        DynType % capacitance coefficient matrix
        D % capacitance coefficient matrix
        B % input mapping matrix
        
    end
    
    methods
        function obj = GraphModel(varargin)
            if nargin == 1
                obj.graph = varargin{1};
                init(obj);
            end
        end
        
%         function ReorderStates(obj,idxNew)
%             obj.Graph.M = obj.Graph.M(idxNew,:); %this will require an update to obj.Graph.E
%             obj.C_coeff = obj.C_coeff(idxNew,:);
%             obj.x_init  = obj.x_init(idxNew);
%             obj.DynType = obj.DynType(idxNew);
%             %             obj.D       = obj.D(idxNew);
%             
%             E = obj.Graph.E;
%             E(:,1) = (1:size(obj.Graph.M,1))*(obj.Graph.M == 1); % set edge matrix tails
%             E(:,2) = (1:size(obj.Graph.M,1))*(obj.Graph.M == -1); % set edge matrix heads
%             obj.Graph.E = E;
%         end

        function init(obj)
            % make vertex matrices
            obj.x_init  = vertcat(obj.graph.Vertices.Initial); 
            obj.DynType = vertcat(obj.graph.Vertices.Type); 
            
            % D matrix
            Dmat = zeros(obj.graph.v_tot,obj.graph.Nee);
            E_idx = arrayfun(@(x) find(x==obj.graph.InternalVertices),vertcat(obj.graph.ExternalEdges.AffectedVertex));            
            for i  = 1:length(E_idx)
                Dmat(E_idx(i),i) = 1;
            end
            obj.D = Dmat;
            
            % B matrix
            Eint = obj.graph.InternalEdges;
            numU = max(arrayfun(@(x) length(x.Input),Eint));
            for j = 1:numU
               obj.B.(['B',num2str(j)]) = zeros(obj.graph.Ne,obj.graph.Nu);
               for i = 1:length(Eint)
                   try
                       obj.B.(['B',num2str(j)])(i,Eint(i).Input(j)) = 1;
                   end
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
            
            
            
            
            
        end
        
        function Modify(obj)
        
        end
        
        function Simulate(obj)
            
        end
        
        
        function SolveGraph(obj) % this function will only work for symbolic expressions at the moment
        
            P = CalcP(Sys,x_full,u); % calculates power flows
            C = CalcC(Sys,x_full); % calcualtes capacitance
            
            eqnA(1:sum(idx_x_a),1) = -Sys.graph.M(idx_x_a,:)*P + Sys.D(idx_x_a,:)*P_e == 0; % system of algebraic equations
            
            [A,Bu] = equationsToMatrix(eqnA,x_a); % convert eqnA to the form Ax=B
            x_a_solution = linsolve(A,Bu); % find solution to the algebraic system
            
            eqnD(1:sum(idx_x_d),1) = diag(C(idx_x_d))^-1*(-Sys.graph.M(idx_x_d,:)*P + Sys.D(idx_x_d,:)*P_e); % system of dynamic equations (
            x_d_solution = subs(eqnD,x_a,x_a_solution); % plug in the algebraic system solution into the dynamic system equations
           
            obj.f = x_d_solution;
            obj.g = [x_full(idx_x_d);x_a_solution]; %y = [x_d; x_a(x_d,x_e,u,P_e)]

            
        end
        
    end
    
    
    
    methods(Static)
        function g_sys = Combine(G, ConnectV, ConnectE) % Create a new GraphModel Object
            % Algorithm 1
            % Algorithm 2
        end
        
    end
        
end

