classdef LookupFunction
    %LOOKUPFUNCTION Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Function (1,1) Type
        Breakpoints (:,1) GraphVertex
    end
    
    methods
        function obj = LookupFunction(varargin)
            if nargin > 1
                obj = my_inputparser(obj,varargin{:});
            end
        end
    end
end

