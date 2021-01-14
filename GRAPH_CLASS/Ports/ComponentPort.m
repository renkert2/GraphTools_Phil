classdef ComponentPort < handle
    %COMPONENTPORT Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Description string
        Element = GraphEdge.empty() %mustBeVertexOrEdge
    end
    
    properties (SetAccess = protected)
        Type PortTypes = 1 % Interconnection Types
        Domain Domains = 'Abstract' % Domain of Port
    end
    
    methods
        function obj = ComponentPort(varargin)
            if nargin == 1
                if isstring(varargin{1}) || ischar(varargin{1})
                    obj.Description = string(varargin{1});
                elseif isa(varargin{1}, 'GraphEdge') || isa(varargin{1}, 'GraphVertex')
                    obj.Element = varargin{1};
                else
                    error('Invalid argument');
                end
            elseif nargin >= 2
                my_inputparser(obj,varargin{:});
            end
        end
        
        function set.Element(obj, val)
            edge_flag = isa(val, 'GraphEdge');
            vertex_flag = isa(val, 'GraphVertex');
            if  edge_flag
                obj.Type = 1;
                verts = [val.TailVertex, val.HeadVertex];
                state_vert_flags = arrayfun(@(x) isa(x,'GraphVertex_Internal'),verts);
                assert(sum(state_vert_flags) == 1, 'Edge Port must connect between an Internal Vertex and an External Vertex');
                state_vert = verts(state_vert_flags);
                obj.Domain = state_vert.VertexType.Domain;
            elseif vertex_flag
                obj.Type = 2;
                if numel(val) == 1
                    obj.Domain = val.VertexType.Domain;
                else
                    state_vert_flags = arrayfun(@(x) isa(x,'GraphVertex_Internal'),val);
                    assert(sum(state_vert_flags) <= 1, 'Multi-vertex Port must contain at most one internal vertex');
                    types = [val.VertexType];
                    assert(all(types(1) == types(2:end)), 'Multi-vertex Port must be of the same VertexType');
                    obj.Domain = types(1).Domain;
                end
            else
                error('Port Element must be a GraphEdge or a GraphVertex')
            end
            obj.Element = val;
        end
    end
end


