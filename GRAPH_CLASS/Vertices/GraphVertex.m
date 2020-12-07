classdef GraphVertex < matlab.mixin.Heterogeneous & matlab.mixin.Copyable
    %GRAPHVERTEX Super Class for all vertex types (i.e. State,
    %External)
    %   The matlab.mixin.Heterogeneous type will allow us to create heterogeneous object
    %   arrays of all the child classes.  
    properties
        Description string = string.empty()
    end
    
    methods
        function obj = GraphVertex(varargin)
            if nargin > 1
                obj = my_inputparser(obj,varargin{:});
            end
        end
    end
end

