classdef GraphEdge < matlab.mixin.Heterogeneous & matlab.mixin.Copyable
    % GraphEdge Super Class for all edge types (i.e. Internal, External) in
    % the Graph Modeling Toolbox
    % Instatiate an empty object or use an input parser
    % Definable properties (with class) include:
    % - Description (string)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Contributors: Christopher T. Aksland and Phil Renkert
    % Association: University of Illionis at Urbana-Champaign
    % Contact: aksland2@illinois.edu and renkert2@illinois.edu
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % potential improvements:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    properties
        Description (1,1) string % Edge Description
    end
    
    properties (SetAccess = ?SystemElement) % Parent can only be set by the instantiating Component
        Parent SystemElement = SystemElement.empty() % edge parent object
    end

    methods
        function obj = GraphEdge(varargin) % constructor method using input parsing
            if nargin > 1
                obj = my_inputparser(obj,varargin{:});
            end
        end
    end
    
    methods (Sealed) % Sealed attribute required for heterogenous arrays
        function x = eq(obj1,obj2) 
            % custom GraphEdge object eq() method to compare a hetergeneous
            % array of GraphEdges
            % GraphEdges are 'handle' objects, so they are considered 
            % equivalent if the handles point to the same GraphEdge object. 
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
            % Custom ne() method; refer to eq() for details
            x = ~eq(obj1, obj2);
        end
        
        function [e,i] = getInternal(obj_array)
            i = arrayfun(@(x) isa(x,'GraphEdge_Internal'),obj_array);
            e = obj_array(i);
        end
        
        function [e,i] = getExternal(obj_array)
            [~,i] = getInternal(obj_array);
            i = ~i;
            e = obj_array(i);
        end
        
        function [t,S] = table(obj_array)
            for i = 1:numel(obj_array)
                v = obj_array(i);
                s = struct();
                s.Number = i;
                s.Parent = v.Parent.Name;
                s.Description = v.Description;
                if isa(v, "GraphEdge_Internal")
                    coeff = v.Coefficient;
                    if numel(coeff) > 1
                        coeff = coeff(1);
                    end
                    flow = v.PowerFlow.Val_Sym;
                    s.PowerFlow = string(coeff*flow);
                    if ~isempty(v.Input)
                        s.Input = v.Input.Description;
                    else
                        s.Input = "---";
                    end
                elseif isa(v, "GraphEdge_External")
                    s.PowerFlow = "---";
                    s.Input = "---";
                else
                    error("Invalid vertex class")
                end
                
                S(i) = s;
            end
            t = struct2table(S);
            t.Properties.VariableNames = {'Number', 'Parent', 'Description', 'Power Flow', 'Input'};
        end
    end
    
end