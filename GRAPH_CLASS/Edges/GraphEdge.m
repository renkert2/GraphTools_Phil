classdef GraphEdge < matlab.mixin.Heterogeneous & matlab.mixin.Copyable
    %GRAPHEdge Super Class for all edge types (i.e. Internal,
    %External)
    %   The matlab.mixin.Heterogeneous type will allow us to create heterogeneous object
    %   arrays of all the child classes.
    
    properties
        Description (1,1) string
        HeadVertex (1,1) GraphVertex = GraphVertex()
    end

    methods
        function obj = GraphEdge(varargin)
            if nargin > 1
                obj = my_inputparser(obj,varargin{:});
            end
        end
        
                
        function set.HeadVertex(obj, head_vertex)
            assert(numel(head_vertex)<=1, 'Edge can have only one head vertex');
            obj.HeadVertex = head_vertex;
        end
    end
    
end