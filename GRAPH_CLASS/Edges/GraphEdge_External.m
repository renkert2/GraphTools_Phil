classdef GraphEdge_External < GraphEdge
    %GRAPHEDGE_External Summary of this class goes here
    %   Detailed explanation goes here
    properties
        HeadVertex GraphVertex = GraphVertex.empty()
    end
    
    methods
        function obj = GraphEdge_External(varargin) % constructor method
            obj@GraphEdge(varargin{:}); % calls the superclass constructor
        end
               
        function set.HeadVertex(obj, head_vertex)
            assert(numel(head_vertex)<=1, 'Edge can have only one head vertex');
            obj.HeadVertex = head_vertex;
        end
    end
    
end

