classdef ComponentPort < handle
    %COMPONENTPORT Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Description string
        Domain DomainTypes = "Abstract"
        Element = GraphEdge.empty() %mustBeVertexOrEdge
    end
    
    properties (SetAccess = protected)
        Type PortTypes = 1 % Interconnection Types
    end
    
    methods
        function obj = ComponentPort(varargin)
            if nargin == 1
                obj.Domain = varargin{1};
            elseif nargin >= 2
                my_inputparser(obj,varargin{:});
            end
        end
        
        function set.Element(obj, val)
            edge_flag = isa(val, 'GraphEdge');
            vertex_flag = isa(val, 'GraphVertex');
            if  edge_flag
                obj.Type = 1;
            elseif vertex_flag
                obj.Type = 2;
            else
                error('Port Element must be a GraphEdge or a GraphVertex')
            end
            
            obj.Element = val;
        end
    end
end


