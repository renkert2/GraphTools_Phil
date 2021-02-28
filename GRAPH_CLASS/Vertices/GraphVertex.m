classdef GraphVertex < handle & matlab.mixin.Heterogeneous & matlab.mixin.Copyable
    % GraphVertex Super Class for all vertex types (i.e. Internal, 
    % External) in the Graph Modeling Toolbox.
    % Instatiate an empty object or use an input parser
    % Definable properties (with class) include:
    % - Description (string)
    % - VertexType (VertexTypes)
    % - DynamicType (DynamicTypes)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Contributors: Christopher T. Aksland and Phil Renkert
    % Association: University of Illionis at Urbana-Champaign
    % Contact: aksland2@illinois.edu and renkert2@illinois.edu
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % potential improvements:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    properties
        Description string = string.empty() % vertex state description
        VertexType VertexTypes = "Abstract" % describes vertex type
        DynamicType DynamicTypes = "EnergyFlow" % Dynamic Type: 1 - Energy Flow, 2 - State Flow
    end
    
    properties (SetAccess = ?Component)
        Parent Component = Component.empty() % parent object for vertices
    end
    
    methods
        function obj = GraphVertex(varargin) % Object constructor
            if nargin > 1
                obj = my_inputparser(obj,varargin{:});
            end
        end
    end
    
    methods (Sealed)   
        function x = eq(obj1,obj2)
            % custom GraphVertex object eq() method to compare a hetergeneous
            % array of GraphVertices
            % GraphVertices are 'handle' objects, so they are considered
            % equivalent if the handles point to the same GraphVertex object.
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
        
        function [v,i] = getInternal(obj_array)
            i = arrayfun(@(x) isa(x,'GraphVertex_Internal'), obj_array);
            v = obj_array(i);
        end
        
        function [v,i] = getExternal(obj_array)
            [~,i] = getInternal(obj_array);
            i = ~i;
            v = obj_array(i);
        end
            
        function [v,i] = getDynamic(obj_array)
            i = arrayfun(@isDynamic,obj_array);
            v = obj_array(i);
            
            function l = isDynamic(x)
                if isa(x,"GraphVertex_Internal")
                    coeff = x.Coefficient;

                    if isa(coeff, 'sym')
                        l = ~ all(arrayfun(@(x) isequal(x, sym(0)), coeff));
                    elseif isa(coeff, 'double')
                        l = sum(abs(coeff))>0;
                    else
                        error('Invalid coefficient type')
                    end
                else
                    l = false;
                end
            end
        end
        
        function [v,i] = getAlgebraic(obj_array)
            [~,i_dyn] = getDynamic(obj_array);
            [~,i_int] = getInternal(obj_array);
            
            i = (~i_dyn) & i_int;
            v = obj_array(i);
        end
                   
    end
end

