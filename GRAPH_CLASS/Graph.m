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
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    properties
        Vertices (:,1) GraphVertex = GraphVertex.empty()
        Edges (:,1) GraphEdge = GraphEdge.empty()
        E % Graph Edge Matrix
    end
    
    properties (SetAccess = private)
        % Incidence Matrix
        M 
        
        % Number of vertices
        Nv
        
        % Number of edges
        Ne
        
        % Number of inputs
        Nu
        
        % Number of external vertices
        Nev
        
        % Number of external edges
        Nee
      
    end
    
    methods
        
        function obj = Graph(varargin) % constructor method
            % this function can be initialized with an edge matrix.
            % Immediately calculate the incidence matrix
           if nargin > 0
               obj.E = varargin{1};
               obj.Vertices = varargin{2};
               obj.Edges = varargin{3};
               
               obj.Nv  = sum(arrayfun(@(x) isa(x,'GraphVertex_Internal'),obj.Vertices));
               obj.Ne  = sum(arrayfun(@(x) isa(x,'GraphEdge_Internal'),obj.Edges));
%                obj.Nu  = ;
               obj.Nev = sum(arrayfun(@(x) isa(x,'GraphVertex_External'),obj.Vertices));
               obj.Nee = sum(arrayfun(@(x) isa(x,'GraphEdge_External'),obj.Edges));
               
               obj.M = zeros(max(max(obj.E)),size(obj.E,1));
               for i = 1:size(obj.E,1)
                   obj.M(obj.E(i,1),i) =  1; % tails
                   obj.M(obj.E(i,2),i) = -1; % heads
               end
           end
               
        end
     
        
        function init(obj) 
            % placeholder
        end
    end
end

