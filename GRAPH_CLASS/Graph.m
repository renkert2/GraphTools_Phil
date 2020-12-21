classdef Graph < matlab.mixin.Copyable
    %GRAPH Fundamental Graph Class
    %   Keep pure and simple: A graph consists of Vertices, Edges, and the
    %   Connection Matrix / Incidence Matrix
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Author: Christopher T. Aksland and Phil Renkert
    % Association: University of Illionis at Urbana-Champaign
    % Contact: aksland2@illinois.edu
    % Revision History:
    % 12/5/2020 - Consider added Matlab "digraph" or adjcancy matrix as a 
    %           propertey.
    % 12/20/2020 - A graph is now entirely defined by edges and vertices.
    %           The edge and incidence matrix are dependent properties.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % potential improvements:
    % Consider added Matlab "digraph" or adjcancy matrix as a propertey.
    % Consider adding dynamic and algebric index get functions
    % Objects being passed in arbitrary order
    % Consider improving how we find the head and tail vertex indices.
    % At some point, graph can be defined with only edges
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    
    properties
        % set of vertices
        Vertices (:,1) GraphVertex = GraphVertex.empty() 
        % set of edges
        Edges (:,1) GraphEdge = GraphEdge.empty() 
    end
    
    properties %(SetAccess = immutable)
            end
    
    properties (SetAccess = private)
    
        % Number of vertices
        Nv (1,1) double %{mustBeInteger,mustBeNonnegative}
        
        % Number of edges
        Ne (1,1) double %{mustBeInteger,mustBeNonnegative}
        
        % Number of inputs
        Nu (1,1) double %{mustBeInteger,mustBeNonnegative}
        
        % Number of external vertices
        Nev (1,1) double %{mustBeInteger,mustBeNonnegative}
        
        % Number of external edges
        Nee (1,1) double %{mustBeInteger,mustBeNonnegative}
        
        % Graph Edge Matrix
        E (:,2) double %{mustBeInteger,mustBePositive}; requires default value for Validation Functions
        
        % Incidence Matrix
        M (:,:) double %{mustBeInteger,mustBeMember({-1,1,0,[]})}; requires default value for Validation Functions

      
    end
    
    properties (Dependent)
        % total number of vertices
        v_tot (1,1) double %{mustBeInteger,mustBeNonnegative}
        
        % Tails
        Tails (1,1) double 
        
        % Heads
        Heads (1,1) double 
        
        % dynamic and Algebraic vertices
%         DynamicVertices (:,1) double 
%         AlgebraicVertices (:,1) double 
    end
    
    properties (Dependent) % make this hidden at some point
        InternalVertices (:,1) GraphVertex = GraphVertex.empty()
        ExternalVertices (:,1) GraphVertex = GraphVertex.empty()
        InternalEdges    (:,1) GraphEdge   = GraphEdge.empty() 
        ExternalEdges    (:,1) GraphEdge   = GraphEdge.empty() 
    end
    
    methods
        
        function obj = Graph(varargin) % constructor method
            % this function can be initialized with an edge matrix.
            % Immediately calculate the incidence matrix
           if nargin == 0
               % do nothing
           elseif nargin == 2
               obj.Vertices = varargin{1};
               obj.Edges = varargin{2};
               obj.init()
               
           else
              error('A Graph object must be initialized as an empty object or as Graph(Edge Matrix or Incidence Matrix, Vertices, Edges).') 
           end
                              
        end
        
        function m = updateM(obj)
            m = zeros(max(max(obj.E)),size(obj.E,1));
               for i = 1:size(obj.E,1)
                   m(obj.E(i,1),i) =  1; % tails
                   m(obj.E(i,2),i) = -1; % heads
               end
        end
        
        function e = updateE(obj)
            e = zeros(size(obj.M,2),2);
            e(:,1) = (1:size(obj.M,1))*(obj.M == 1); % set edge matrix tails
            e(:,2) = (1:size(obj.M,1))*(obj.M == -1); % set edge matrix heads
        end
            
        function x = get.v_tot(obj)
            x = obj.Nv + obj.Nev;
        end
        
        function x = get.Tails(obj)
            x = double(obj.M'== 1);
        end
        
        function x = get.Heads(obj)
            x = double(obj.M'== -1);
        end
        
        function x = get.InternalVertices(obj) % these should be functions of the vertex and edge objects
            x = obj.Vertices(arrayfun(@(x) isa(x,'GraphVertex_Internal'),obj.Vertices));
        end     
        function x = get.InternalEdges(obj)
            x = obj.Edges(arrayfun(@(x) isa(x,'GraphEdge_Internal'),obj.Edges));
        end
        function x = get.ExternalVertices(obj)
            x = obj.Vertices(~arrayfun(@(x) isa(x,'GraphVertex_Internal'),obj.Vertices));
        end     
        function x = get.ExternalEdges(obj)
            x = obj.Edges(~arrayfun(@(x) isa(x,'GraphEdge_Internal'),obj.Edges));
        end
        
        
        function init(obj) 
            obj.Nv  = length(obj.InternalVertices);
            obj.Ne  = length(obj.InternalEdges);
            obj.Nu  = max([obj.InternalEdges.Input]);
            obj.Nev = length(obj.Vertices)-obj.Nv;
            obj.Nee = length(obj.Edges)-obj.Ne;
            
            vTail = vertcat(obj.InternalEdges.TailVertex);
            vHead = vertcat(obj.InternalEdges.HeadVertex);

            % this finds tails and heads vertices. The code is clunky
            % because we can't use == for a hetergenous class. We could
            % look into speeding this up at some point.
%             vT_idx = arrayfun(@(x) (find(x==obj.Vertices)),vTail);
            
            vT_idx =  arrayfun(@(x) find(arrayfun(@(y) eq(x,y),obj.Vertices)),vTail);
            vH_idx =  arrayfun(@(x) find(arrayfun(@(y) eq(x,y),obj.Vertices)),vHead);

%             isequal(vTail(1),obj.Vertices)
%             vT_idx = arrayfun(@(x) (my_find_eq(x,obj.InternalVertices,0))+(my_find_eq(x,obj.ExternalVertices,obj.Nv)),vTail);
%             vH_idx = arrayfun(@(x) (my_find_eq(x,obj.InternalVertices,0))+(my_find_eq(x,obj.ExternalVertices,obj.Nv)),vHead);

            obj.E = [vT_idx vH_idx];
            m = zeros(max(max(obj.E)),size(obj.E,1));
            for i = 1:size(obj.E,1)
                m(obj.E(i,1),i) =  1; % tails
                m(obj.E(i,2),i) = -1; % heads
            end
            obj.M = m;
            
        end
        
        
        function plot(obj,varargin)
            % basic digraph plotting.
%             figure
            E = obj.E; % edge matrix
            try
                Edge_ext = vertcat(obj.ExternalEdges.HeadVertex);
                E_idx = arrayfun(@(x) find(x==obj.InternalVertices),Edge_ext);
                Eext = [[obj.v_tot+1:1:obj.v_tot+length(E_idx)]' E_idx];
                E = [E; Eext]; % augment E matrix with external edges
                skipPlotExt = 0;
            catch
                skipPlotExt = 1;
            end
            G = digraph(E(:,1),E(:,2));
            h = plot(G,varargin{:});
            labeledge(h,E(:,1)',E(:,2)',[1:obj.Ne, 1:obj.Nee]);
            highlight(h,[obj.Nv+1:1:obj.v_tot],'NodeColor','w')
            xLoc = h.XData(obj.Nv+1:1:obj.v_tot);
            yLoc = h.YData(obj.Nv+1:1:obj.v_tot);
            hold on; scatter(xLoc,yLoc,5*h.MarkerSize,'MarkerEdgeColor',h.NodeColor(1,:)); hold off;
            if ~skipPlotExt
                highlight(h,reshape(Eext',1,[]),'LineStyle','--')
                highlight(h,[obj.v_tot+1:1:obj.v_tot+length(E_idx)],'NodeLabelColor','w')
                highlight(h,[obj.v_tot+1:1:obj.v_tot+length(E_idx)],'NodeColor','w')
                xLoc = [ h.XData(Eext(:,1))];
                yLoc = [ h.YData(Eext(:,1))];
            end
            hold on; scatter(xLoc,yLoc,5*h.MarkerSize,'MarkerEdgeColor',h.EdgeColor(1,:)); hold off;

        end
        
        function ReorderVertices(obj,idx)
            obj.Vertices = obj.Vertices(idx);
            obj.init();
        end
        
    end
end

