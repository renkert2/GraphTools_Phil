classdef ComponentPort < handle
    % The ComponentPort class defines external connections for a component
    % in the Graph Modeling Toolbox. The external connections are
    % associated with an edge of vertex (see PortTypes) of a Graph object.
    % Instatiate an empty object, use the input parser, or
    % P = ComponentPort(desc) defines port P with Description desc or
    % P = ComponentPort(elmt) defines port P with Element elmt
    % Definable properties (with class) include:
    % - Description (string)
    % - Element (GraphEdge or GraphVertex)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Contributors: Christopher T. Aksland and Phil Renkert
    % Association: University of Illionis at Urbana-Champaign
    % Contact: aksland2@illinois.edu and renkert2@illinois.edu
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % potential improvements:
    %  - remove the else state block in set.Element
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    properties
        Description (1,1) string
        Element (1,1) %= GraphEdge.empty() %mustBeVertexOrEdge
    end
    
    properties (SetAccess = protected)
        Type (1,1) PortTypes = 1 % Interconnection Types
        Domain (1,1) Domains = 'Abstract' % Domain of Port
    end
    
    methods
        function obj = ComponentPort(varargin) % constructor. Uses input parsing or initializes the object with a description OR element
            if nargin == 1
                if isstring(varargin{1}) || ischar(varargin{1}) % initialize with description
                    obj.Description = string(varargin{1});
                elseif isa(varargin{1}, 'GraphEdge') || isa(varargin{1}, 'GraphVertex') % initialize with element
                    obj.Element = varargin{1};
                else
                    error('Invalid argument'); % throw error
                end
            elseif nargin >= 2 % input parsing
                my_inputparser(obj,varargin{:});
            end
        end
        
        function set.Element(obj, val)
            edge_flag = isa(val, 'GraphEdge');
            vertex_flag = isa(val, 'GraphVertex'); % error check to make sure the element is either a GraphEdge or Graph Vertex object
            if  edge_flag
                obj.Type = 1; % verify that an edge port includes at least 1 GraphVertex_External object
                verts = [val.TailVertex, val.HeadVertex];
                state_vert_flags = arrayfun(@(x) isa(x,'GraphVertex_Internal'),verts);
                assert(sum(state_vert_flags) == 1, 'Edge Port must connect an Internal Vertex and an External Vertex');
                state_vert = verts(state_vert_flags);
                obj.Domain = state_vert.VertexType.Domain;
            elseif vertex_flag
                obj.Type = 2; % verify that an edge port includes at least 1 GraphVertex_External object
                if numel(val) == 1
                    obj.Domain = val.VertexType.Domain;
                else
                    % @ PHIL... A single port can only have a single vertex
                    % or edge. we should not need this else state... Unless
                    % this block is called somewhere else?
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


