classdef GraphEdge_Internal < GraphEdge
    % GraphEdge_Internal defines the properties of internal edges in a
    % graph object in the Graph Modeling Toolbox
    % Instatiate an empty object or use an input parser
    % Definable properties (with class) include:
    % - Description (string)
    % - PowerFlow (Type_PowerFlow) 
    % - Input (GraphInput)
    % - Coefficient (double)
    % - TailVertex (GraphVertex) 
    % - HeadVertex (GraphVertex) 
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Contributors: Christopher T. Aksland and Phil Renkert
    % Association: University of Illionis at Urbana-Champaign
    % Contact: aksland2@illinois.edu and renkert2@illinois.edu
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % potential improvements:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    properties
        PowerFlow (:,1) Type_PowerFlow %= Type_PowerFlow.empty() % powerflow representation for an edge
        Input (:,1) GraphInput %= GraphInput.empty() % inputs incident on the edge 
        Coefficient (:,1) double = 0 % the constant edge coefficient
        TailVertex (1,1) GraphVertex %= GraphVertex.empty() % the edge tail vertex
        HeadVertex (1,1) GraphVertex %= GraphVertex.empty() % the edge head vertex
    end
    
    methods
        function obj = GraphEdge_Internal(varargin) % constructor method
            obj@GraphEdge(varargin{:}); % calls the superclass constructor
        end
        
%         function set.TailVertex(obj, tail_vertex) % error check on whether two vertices can be connected by an edge
%             if ~isempty(obj.HeadVertex)
%                 if ~isConnectable(obj.HeadVertex.VertexType, tail_vertex.VertexType)
%                     warning('Tail vertex is not connectable with existing head vertex')
%                 end
%             end
%             obj.TailVertex = tail_vertex;
%         end
%         
%         function set.HeadVertex(obj, head_vertex) % error check on whether two vertices can be connected by an edge
%             if ~isempty(obj.TailVertex)
%                 if ~isConnectable(obj.TailVertex.VertexType, head_vertex.VertexType)
%                     warning('Head vertex is not connectable with existing tail vertex')
%                 end
%             end
%             obj.HeadVertex = head_vertex;
%         end
    end
end
