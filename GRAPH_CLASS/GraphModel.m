classdef GraphModel < matlab.mixin.Copyable
    %GRAPHMODEL Contains all of the operations useful for working with the
    %graph models.  Most of the code will go here.  
    %   Detailed explanation goes here
    properties
        Graph Graph = Graph.empty()
        
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
                obj.Graph = varargin{1};
            end
        end
        
        function plot(obj,varargin)
            % basic digraph plotting.
%             figure
            E = obj.Graph.E; % edge matrix
            try
                Edge_ext = vertcat(obj.Graph.ExternalEdges.V_ind);
                Edge_ext(Edge_ext == 0) = [];
                Eext = [[obj.Graph.v_tot+1:1:obj.Graph.v_tot+length(Edge_ext)]' Edge_ext];
                E = [E; Eext]; % augment E matrix with external edges
                skipPlotExt = 0;
            catch
                skipPlotExt = 1;
            end
            G = digraph(E(:,1),E(:,2));
            h = plot(G,varargin{:});
            labeledge(h,E(:,1)',E(:,2)',[1:obj.Graph.Ne, 1:obj.Graph.Nee]);
            highlight(h,[obj.Graph.Nv+1:1:obj.Graph.v_tot],'NodeColor','w')
            xLoc = h.XData(obj.Graph.Nv+1:1:obj.Graph.v_tot);
            yLoc = h.YData(obj.Graph.Nv+1:1:obj.Graph.v_tot);
            hold on; scatter(xLoc,yLoc,5*h.MarkerSize,'MarkerEdgeColor',h.NodeColor(1,:)); hold off;
            if ~skipPlotExt
                highlight(h,reshape(Eext',1,[]),'LineStyle','--')
                highlight(h,[obj.Graph.v_tot+1:1:obj.Graph.v_tot+length(Edge_ext)],'NodeLabelColor','w')
                highlight(h,[obj.Graph.v_tot+1:1:obj.Graph.v_tot+length(Edge_ext)],'NodeColor','w')
                xLoc = [ h.XData(Eext(:,1))];
                yLoc = [ h.YData(Eext(:,1))];
            end
            hold on; scatter(xLoc,yLoc,5*h.MarkerSize,'MarkerEdgeColor',h.EdgeColor(1,:)); hold off;

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
            obj.x_init  = vertcat(obj.Graph.Vertices.Initial); 
            obj.DynType = vertcat(obj.Graph.Vertices.Type); 
            
            % D matrix
            Dmat = zeros(obj.Graph.v_tot,obj.Graph.Nee);
            Eext = obj.Graph.ExternalEdges;
            for i  = 1:length(Eext)
                Dmat(Eext(i).V_ind,i) = 1;
            end
            obj.D = Dmat;
            
            % B matrix
            Eint = obj.Graph.InternalEdges;
            numU = max(arrayfun(@(x) length(x.Input),Eint));
            for j = 1:numU
               obj.B.(['B',num2str(j)]) = zeros(obj.Graph.Ne,obj.Graph.Nu);
               for i = 1:length(Eint)
                   try
                       obj.B.(['B',num2str(j)])(i,Eint(i).Input(j)) = 1;
                   end
               end
            end
            
            % C matrix
            CTypeAll = vertcat(obj.Graph.Vertices(:).Capacitance);
            numCType = arrayfun(@(x) length(x.Capacitance),obj.Graph.Vertices);
            [obj.C_coeff,obj.CType] = MakeCoeffMatrix(obj.Graph.Vertices,CTypeAll,numCType);

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

