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
    % 12/21/2020 - A graph edge and vertex set can now be defined in
    %           arbitrary order.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % potential improvements:
    % Consider added Matlab "digraph" or adjcancy matrix as a propertey.
    % Consider improving how we find the head and tail vertex indices.
    % At some point, graph can be defined with only edges
    % Move the Internal/External etc get functions into the GraphVertex and
    %   GraphEdge classes. Since those functions only operate on vertices
    %   and edges, it seems they fit better there.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    
    properties
        % set of vertices
        Vertices (:,1) GraphVertex = GraphVertex.empty() 
        % set of edges
        Edges (:,1) GraphEdge = GraphEdge.empty() 
    end
    
    properties (SetAccess = private)
        Inputs (:,1) GraphInput = GraphInput.empty()
    
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
        
    end
    
    properties (Dependent) % make this hidden at some point
        InternalVertices  (:,1) GraphVertex = GraphVertex.empty()
        DynamicVertices   (:,1) GraphVertex = GraphVertex.empty()
        AlgebraicVertices (:,1) GraphVertex = GraphVertex.empty()
        ExternalVertices  (:,1) GraphVertex = GraphVertex.empty()
        InternalEdges     (:,1) GraphEdge   = GraphEdge.empty() 
        ExternalEdges     (:,1) GraphEdge   = GraphEdge.empty() 
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
              error('A Graph object must be initialized as an empty object or as Graph(Vertex_Set, Edge_Set ).') 
           end
                              
        end
                      
        function x = get.v_tot(obj) % total number of vertices
            x = obj.Nv + obj.Nev;
        end
        
        function x = get.Tails(obj) % tail incidence matrix
            x = double(obj.M'== 1);
        end
        
        function x = get.Heads(obj) % head incidence matrix
            x = double(obj.M'== -1);
        end
        
        function x = get.InternalVertices(obj) % these should be functions of the vertex and edge objects?
            x = obj.Vertices(arrayfun(@(x) isa(x,'GraphVertex_Internal'),obj.Vertices));
        end 
        function x = get.DynamicVertices(obj) % these should be functions of the vertex and edge objects
            x = obj.Vertices(arrayfun(@(x) (isa(x,'GraphVertex_Internal') && (sum(abs(x.Coefficient))>0)),obj.Vertices));
        end
        function x = get.AlgebraicVertices(obj) % these should be functions of the vertex and edge objects
            x = obj.Vertices(arrayfun(@(x) (isa(x,'GraphVertex_Internal') && (sum(abs(x.Coefficient))==0)),obj.Vertices));
        end
        function x = get.InternalEdges(obj) % these should be functions of the vertex and edge objects
            x = obj.Edges(arrayfun(@(x) isa(x,'GraphEdge_Internal'),obj.Edges));
        end
        function x = get.ExternalVertices(obj) % these should be functions of the vertex and edge objects
            x = obj.Vertices(~arrayfun(@(x) isa(x,'GraphVertex_Internal'),obj.Vertices));
        end     
        function x = get.ExternalEdges(obj) % these should be functions of the vertex and edge objects
            x = obj.Edges(~arrayfun(@(x) isa(x,'GraphEdge_Internal'),obj.Edges));
        end
        
        
        function init(obj)
            
            % reorder vertices in order of [xd;xa;xt] and edges in [int; ext]
            obj.Vertices = [obj.DynamicVertices; obj.AlgebraicVertices; obj.ExternalVertices];
            obj.Edges    = [obj.InternalEdges; obj.ExternalEdges];            
            
            % Calculate input array from Internal Edges
            obj.Inputs = unique(vertcat(obj.InternalEdges.Input), 'stable');

            % calculate graph size, order, etc
            obj.Nv  = length(obj.InternalVertices);
            obj.Ne  = length(obj.InternalEdges);
            obj.Nu  = length(obj.Inputs);
            obj.Nev = length(obj.Vertices)-obj.Nv;
            obj.Nee = length(obj.Edges)-obj.Ne;
            
            % list of tail and head vertices
            vTail = vertcat(obj.InternalEdges.TailVertex);
            vHead = vertcat(obj.InternalEdges.HeadVertex);
            
            % indicies of tail and vertices
            vT_idx =  arrayfun(@(x) find(arrayfun(@(y) eq(x,y),obj.Vertices)),vTail);
            vH_idx =  arrayfun(@(x) find(arrayfun(@(y) eq(x,y),obj.Vertices)),vHead);

            % calculate Edge and Incidence matrix
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
            E_Augmented = obj.E; % edge matrix
            try % try to augment the edge matrix to include external edges
                E_idx = arrayfun(@(x) find(x==obj.InternalVertices),vertcat(obj.ExternalEdges.HeadVertex)); % index of vertices with external edges
                Eext = [[obj.v_tot+1:1:obj.v_tot+obj.Nee]' E_idx]; % edge matrix of external edges
                E_Augmented = [E_Augmented; Eext]; % augment E matrix with external edges
                skipPlotExt = 0;
            catch
                skipPlotExt = 1;
            end
            G = digraph(E_Augmented(:,1),E_Augmented(:,2)); % make digraph
            h = plot(G,varargin{:}); % plot digraph
            labeledge(h,E_Augmented(:,1)',E_Augmented(:,2)',[1:obj.Ne, 1:obj.Nee]); % edge lables
            highlight(h,[obj.Nv+1:1:obj.v_tot],'NodeColor','w') % remove external vertices
            xLoc = h.XData(obj.Nv+1:1:obj.v_tot); yLoc = h.YData(obj.Nv+1:1:obj.v_tot); % get external vertex location
            hold on; scatter(xLoc,yLoc,5*h.MarkerSize,'MarkerEdgeColor',h.NodeColor(1,:)); hold off; % add external vertices back
            if ~skipPlotExt
                highlight(h,reshape(Eext',1,[]),'LineStyle','--') % make external edges dashed
                highlight(h,[obj.v_tot+1:1:obj.v_tot+length(E_idx)],'NodeLabelColor','w','NodeColor','w') % remove augmented external edge tail vertices and labels
                xLoc = [ h.XData(Eext(:,1))]; yLoc = [ h.YData(Eext(:,1))]; % get external edge tail vertex location
            end
            hold on; scatter(xLoc,yLoc,5*h.MarkerSize,'MarkerEdgeColor',h.EdgeColor(1,:)); hold off; % plot external edge tail vertices
        end
    end
end

