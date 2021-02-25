classdef GraphOutput < handle
    %GRAPHOUTPUT Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Description (1,1) string = "Default"
        Function % Must be symfun or Type, asserted in set method
        Breakpoints cell % Must be cell array of GraphVertices or GraphInputs, asserted in set method
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
        
        function set.Function(obj, val)
            assert(isa(val, 'Type') || isa(val, 'symfun'), 'Function must be a Type object or symfun');
            obj.Function = val
        end
        
        function set.Breakpoints(obj,val)
            for i = 1:length(val)
                if ~isa(val{i},'GraphVertex') && ~isa(val{i},'GraphInput')
                    error('Breakpoints must be of type GraphVertex or GraphInput')
                end
            end
            obj.Breakpoints = val;
        end
        
    end
end

