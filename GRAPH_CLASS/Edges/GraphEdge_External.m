classdef GraphEdge_External < GraphEdge
    %GRAPHEDGE_External Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Description (1,1) string
        
    end

    
    
    methods
        function obj = GraphEdge_External(varargin) % constructor method
            obj@GraphEdge(varargin{:}); % calls the superclass constructor
        end       
    end
    
end

