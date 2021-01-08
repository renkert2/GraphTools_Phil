classdef GraphEdge < matlab.mixin.Heterogeneous & matlab.mixin.Copyable
    %GRAPHEdge Super Class for all edge types (i.e. Internal,
    %External)
    %   The matlab.mixin.Heterogeneous type will allow us to create heterogeneous object
    %   arrays of all the child classes.
    
    properties
        Description (1,1) string
        HeadVertex (1,1) GraphVertex = GraphVertex()
    end
    
    properties (SetAccess = ?Component)
        Parent Component = Component.empty()
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
    
    methods (Sealed)
        function x = eq(obj1,obj2)
            if numel(obj1)==numel(obj2)
                x = false(size(obj1));
                for i = 1:numel(obj1)
                    x(i) = eq@handle(obj1(i), obj2(i));
                end
            elseif numel(obj1)==1
                x = false(size(obj2));
                for i = 1:numel(obj2)
                    x(i) = eq@handle(obj1, obj2(i));
                end
            elseif numel(obj2)==1
                x = false(size(obj1));
                for i = 1:numel(obj1)
                    x(i) = eq@handle(obj2, obj1(i));
                end
            else
                error('array sizes are incompatible') 
            end
        end
        
        function x = ne(obj1, obj2)
            x = ~eq(obj1, obj2);
        end
    end
    
end