classdef GraphInput < handle
    %GraphInput
    
    properties
        Description string = string.empty()
    end
    
    methods
        function obj = GraphInput(varargin)
            if nargin > 1
                obj = my_inputparser(obj,varargin{:});
            end
        end
    end
end

