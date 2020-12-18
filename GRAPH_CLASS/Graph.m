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
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % potential improvements:
    % Consider added Matlab "digraph" or adjcancy matrix as a propertey.
    % consider adding dynamic and algebric index get functions
    % note that M and E are technically dependent on Vertices and Edges...
    % at some point, we'll need to account for internal and external
    % objects being passed in arbitrary order
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    
    properties
        % set of vertices
        Vertices (:,1) GraphVertex = GraphVertex.empty() 
        % set of edges
        Edges (:,1) GraphEdge = GraphEdge.empty() 
    end
    
    properties %(SetAccess = immutable)
        % Graph Edge Matrix
        E (:,2) double %{mustBeInteger,mustBePositive}; requires default value for Validation Functions
        % Incidence Matrix
        M (:,:) double %{mustBeInteger,mustBeMember({-1,1,0,[]})}; requires default value for Validation Functions
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
           elseif nargin == 3
               if any(varargin{1}<0)
                   obj.M = varargin{1}; % define Incidence matrix
                   obj.E = updateE(obj);
               else
                   obj.E = varargin{1}; % define Edge Matrix
                   obj.M = updateM(obj); % define Edge Matrix
               end
               obj.Vertices = varargin{2};
               obj.Edges = varargin{3};
               
               obj.Nv  = length(obj.InternalVertices);
               obj.Ne  = length(obj.InternalEdges);
               obj.Nu  = max([obj.InternalEdges.Input]);
               obj.Nev = length(obj.Vertices)-obj.Nv;
               obj.Nee = length(obj.Edges)-obj.Ne;
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
        
        function x = get.InternalVertices(obj)
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
            % placeholder
        end
    end
end

