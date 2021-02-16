classdef GraphOutput < handle
    %GRAPHOUTPUT Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Description (1,1) string = "Default"
        Function (1,1) Type
        Breakpoints 
    end
    
    properties (SetAccess = ?Component)
        Parent Component = Component.empty()
    end
    
    methods
        function obj = GraphOutput(varargin)
            if nargin > 1
                obj = my_inputparser(obj,varargin{:});
            end
        end
        
        function obj = set.Breakpoints(obj,val)
            for i = 1:length(val)
                if ~isa(val{i},'GraphVertex') && ~isa(val{i},'GraphInput')
                    error('Breakpoints must be of type GraphVertex or GraphInput')
                end
            end
            obj.Breakpoints = val;
        end
        
    end
end

