classdef GraphVertex_Internal < GraphVertex
    %GRAPHVERTEX_INTERNAL Summary of this class goes here
    %   Detailed explanation goes here  
    properties
        Type % Vertex Type: 1 - Energy Flow, 2 - State Flow
        Capacitance %Type_Capacitance
    end
    
    methods
        function obj = GraphVertex_Internal(varargin) % constructor method
            obj@GraphVertex(varargin{:}); % calls the superclass constructor
        end       
    end
end

