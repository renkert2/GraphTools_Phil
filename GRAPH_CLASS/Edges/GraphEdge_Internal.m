classdef GraphEdge_Internal < GraphEdge
    %GraphEdge_Internal Summary of this class goes here
    %   Detailed explanation goes here
    properties
        PowerFlow Type_PowerFlow = Type_PowerFlow.empty()
        Input
        Port
    end
    
    methods
        function obj = GraphEdge_Internal(varargin) % constructor method
            obj@GraphEdge(varargin{:}); % calls the superclass constructor
        end       
    end
end

