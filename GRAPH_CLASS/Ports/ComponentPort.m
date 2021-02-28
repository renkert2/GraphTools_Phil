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
        Element (:,1) =  GraphEdge.empty() %must be Vertex Or Edge, asserted in set method
    end
    
    properties (SetAccess = protected)
        Type (1,1) PortTypes = 1 % Interconnection Types: Type 1 for edge connections and Type 2 for vertex connections
        Domain (1,1) Domains = 'Abstract' % Domain of Port
    end
    
    methods
        function obj = ComponentPort(varargin) % Constructor. Uses input parsing or initializes the object with a description OR element
            if nargin == 1
                if isstring(varargin{1}) || ischar(varargin{1}) % initialize with description
                    obj.Description = string(varargin{1});
                elseif isa(varargin{1}, 'GraphEdge') || isa(varargin{1}, 'GraphVertex') % initialize with element
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
                assert(sum(state_vert_flags) == 1, 'Edge Port must connect an Internal Vertex and an External Vertex');
                state_vert = verts(state_vert_flags);
                obj.Domain = state_vert.VertexType.Domain;
            elseif vertex_flag
                obj.Type = 2;
                if numel(val) == 1
                    % Single Vertex Port
                    
                    obj.Domain = val.VertexType.Domain;
                else
                    % Multi-Vertex Port: Used to combine multiple graph vertices.  
                    % - All elements must be connectable.
                    % - Vertices in multi-vertex port will be combined into single vertex after calling Component.Combine()
                    
                    warning("Defining a Multi-Vertex Port.  All elements of a multi-vertex port will be combined into a single vertex when the component is combined with others in Component.Combine()")
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


