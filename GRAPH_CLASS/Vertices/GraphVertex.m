classdef GraphVertex < matlab.mixin.Heterogeneous & matlab.mixin.Copyable
    %GRAPHVERTEX Super Class for all vertex types (i.e. State,
    %External)
    %   The matlab.mixin.Heterogeneous type will allow us to create heterogeneous object
    %   arrays of all the child classes.  
    properties
        Description string = string.empty()
        Type = 1 % Vertex Type: 1 - Energy Flow, 2 - State Flow
        Capacitance Type_Capacitance = Type_Capacitance.empty();
        Coefficient = 0
        Initial = 0
    end
    
    methods
        function obj = GraphVertex(varargin)
            if nargin > 1
                obj = my_inputparser(obj,varargin{:});
            end
        end
    end
end

