classdef GraphVertex_Internal < GraphVertex
    %GRAPHVERTEX_INTERNAL Summary of this class goes here
    %   Detailed explanation goes here  
    properties
        CapFunction (1,1) LookupFunction
        Capacitance (:,1) Type_Capacitance = Type_Capacitance.empty();
        Coefficient (:,1) {mustBeNonnegative} = 0
        Initial (1,1) double = 0
        Bounds Limits = Limits()
    end
    
    methods
        function obj = GraphVertex_Internal(varargin) % constructor method
            obj@GraphVertex(varargin{:}); % calls the superclass constructor
        end       
    end
end

