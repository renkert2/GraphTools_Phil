classdef ComponentPort < handle
    %COMPONENTPORT Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Description string
        PortType PortTypes
        Edge GraphEdge
    end
    
    methods
        function obj = ComponentPort(varargin)
            if nargin == 1
                obj.PortType = varargin{1};
            elseif nargin >= 2
                my_inputparser(obj,varargin{:});
            end
        end
    end
end

