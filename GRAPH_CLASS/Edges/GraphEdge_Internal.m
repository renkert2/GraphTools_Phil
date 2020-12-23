classdef GraphEdge_Internal < GraphEdge
    %GraphEdge_Internal Summary of this class goes here
    %   Detailed explanation goes here
    properties
        PowerFlow Type_PowerFlow = Type_PowerFlow.empty()
        Input GraphInput = GraphInput.empty()
        Port
        Coefficient double = 0
        TailVertex GraphVertex = GraphVertex()

    end
    
    methods
        function obj = GraphEdge_Internal(varargin) % constructor method
            obj@GraphEdge(varargin{:}); % calls the superclass constructor
        end
        
        function set.TailVertex(obj, tail_vertex)
            assert(numel(tail_vertex)<=1, 'Edge can have only one tail vertex');
            obj.TailVertex = tail_vertex;
        end
    end
end
