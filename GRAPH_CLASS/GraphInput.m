classdef GraphInput < handle
    %GraphInput
    
    properties
        Description string = string.empty() 
        Bounds Limits = Limits() 
    end
    
    properties (SetAccess = ?Component)
        Parent Component = Component.empty()
    end
    
    methods
        function obj = GraphInput(varargin)
            if nargin == 1 && (isstring(varargin{1}) || ischar(varargin{1}))
                desc = string(varargin{1});
                obj.Description = desc;
            elseif nargin > 1
                obj = my_inputparser(obj,varargin{:});
            else
                error('Invalid argument to GraphInput constructor');
            end
        end
    end
end

