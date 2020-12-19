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
                MakeMatrices(obj);
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

        function MakeMatrices(obj)
            % make vertex matrices
            obj.x_init  = vertcat(obj.graph.Vertices.Initial); 
            obj.DynType = vertcat(obj.graph.Vertices.Type); 
            
            % D matrix
            Dmat = zeros(obj.graph.v_tot,obj.graph.Nee);
            Eext = obj.graph.ExternalEdges;
            for i  = 1:length(Eext)
                Dmat(Eext(i).V_ind,i) = 1;
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
        
        function init(obj)
            % placeholder
        end
    end
    
    
    
    methods(Static)
        function g_sys = Combine(G, ConnectV, ConnectE) % Create a new GraphModel Object
            % Algorithm 1
            % Algorithm 2
        end
        
    end
        
end

