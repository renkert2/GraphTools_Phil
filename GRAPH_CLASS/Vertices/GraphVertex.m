classdef GraphVertex < handle & matlab.mixin.Heterogeneous & matlab.mixin.Copyable
    %GRAPHVERTEX Super Class for all vertex types (i.e. State,
    %External)
    %   The matlab.mixin.Heterogeneous type will allow us to create heterogeneous object
    %   arrays of all the child classes.  
    properties
        Description string = string.empty()
        VertexType VertexTypes = "Abstract"
        DynamicType DynamicTypes = "EnergyFlow" % Dynamic Type: 1 - Energy Flow, 2 - State Flow
        Capacitance (:,1) Type_Capacitance = Type_Capacitance.empty();
        Coefficient (:,1) = 0
        Initial = 0
    end
    
    properties (SetAccess = ?Component)
        Parent Component = Component.empty()
    end
    
    methods
        function obj = GraphVertex(varargin)
            if nargin > 1
                obj = my_inputparser(obj,varargin{:});
            end
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

