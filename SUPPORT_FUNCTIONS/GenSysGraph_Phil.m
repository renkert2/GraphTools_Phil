function [G_sys, ConnectE, ConnectV] = GenSysGraph_Phil(G, ConnectE, ConnectV, opts) % Input is vector of GraphClass elements
% G - GraphClass Class Array
% ConnectX - Cell vector with first list containing component
% indices; second list containing equivalent properties (edges
% or vertices) corresponding to the components.

arguments
    G (:,1) Graph
    ConnectE cell = {}
    ConnectV cell = {}
    opts.CopyEdges logical = true;
end

% Format ConnectV and ConnectE as {:,1} cell array containing list of
% equivalent vertices and edges, respectively
ConnectE = formatConnectX(ConnectE, 'GraphEdge', 'Edges');
Nce = length(ConnectE); % Number of edge connections
Neconn_delta = cellfun(@length, ConnectE);% Number of edges involved in each connection
Neconn = sum(Neconn_delta); % Total edges involved in connections

if nargin == 3 % If ConnectV is specified
    ConnectV = formatConnectX(ConnectV, 'GraphVertex', 'Vertices');
end

%% Parse ConnectE for necessary Vertex Connections and create edge_conn_map:
% - Appends ConnectV with necessary vertex connections resulting from edge connections
% edge_conn_map:
% - First column contains edges replaced in connection that will be modified in the connection, 
% - Second column contains primary edges resulting from connection
% - s.t. edge_conn_map(:,1) becomes edge_conn_map(:,2)
ConnectV_E = cell(2*Nce, 1);
edge_conn_map(Nce,2) = GraphEdge();
for ce = 1:Nce
    edges = fliplr(ConnectE{ce});

    equiv_heads = [edges.HeadVertex];
    equiv_tails = [edges.TailVertex];
    
    ConnectV_E{2*ce - 1} = equiv_heads;
    ConnectV_E{2*ce} = equiv_tails;
    
    edge_conn_map(ce, :) = edges;
end

if isempty(ConnectV)
    ConnectV = ConnectV_E;
else
    ConnectV = [ConnectV; ConnectV_E];
end
%% Create vertex_conn_map:
% - First column contains all vertices that will be modified in the connection, 
% - Second column contains primary vertices that verts in the first column will map to
% - s.t. Vertex_Conn_Map(:,1) becomes Vertex_Conn_Map(:,2)

% Determine the Primary Vertex for each connection
Ncv = length(ConnectV); % Number of vertex connections
Nvconn_delta = cellfun(@length, ConnectV); % Number of vertices involved in each connection
Nvconn = sum(Nvconn_delta); % Total vertices involved in connections

primary_vertices(Ncv,1) = GraphVertex();
vertex_conn_map(Nvconn-Ncv, 2) = GraphVertex();
vconnmap_counter = cumsum([0; Nvconn_delta-1]);
for cv = 1:length(ConnectV)
    verts = ConnectV{cv};
    int_vert_flags = arrayfun(@(x) isa(x,'GraphVertex_Internal'), verts);
    n_int_verts = sum(int_vert_flags);
    if n_int_verts == 0
         primary_vertex = verts(1);
    elseif n_int_verts == 1
        primary_vertex = verts(int_vert_flags);
    else
        error('Error in Vertex Connection %d.  Only one Internal Vertex can be involved in a connection.', cv)
    end
    primary_vertices(cv,1) = primary_vertex;
    
    conn_verts = verts(primary_vertex ~= verts);
    
    r = vconnmap_counter(cv)+(1:numel(conn_verts)); % Range of elements in vertex_conn_map to be assigned
    vertex_conn_map(r, 1) = conn_verts;
    vertex_conn_map(r,2) = primary_vertex;
end

vertex_conn_map = unique(vertex_conn_map, 'rows','stable'); % remove duplicate mappings in vertex_conn_map
if length(vertex_conn_map(:,1)) ~= length(unique(vertex_conn_map(:,1)))
    error('Vertices assigned to multiple primary vertices.');
end

%% Construct Vertex and Edge Vectors
all_verts = vertcat(G.Vertices);
all_edges = vertcat(G.Edges);

sys_verts = all_verts(~ismember(all_verts, vertex_conn_map(:,1)));
sys_edges = all_edges(~ismember(all_edges, edge_conn_map(:,1)));

if opts.CopyEdges
    sys_edges = copy(sys_edges);
end

for i = 1:numel(sys_edges)
    head_v = sys_edges(i).HeadVertex;
    
    log_index = head_v == vertex_conn_map(:,1);
    if any(log_index)
        primary_vertex = vertex_conn_map(log_index, 2);
        sys_edges(i).HeadVertex = primary_vertex;
    end
    
    if isa(sys_edges(i), 'GraphEdge_Internal')
        tail_v = sys_edges(i).TailVertex;
        
        log_index = tail_v == vertex_conn_map(:,1);
        if any(log_index)
            primary_vertex = vertex_conn_map(log_index, 2);
            sys_edges(i).TailVertex = primary_vertex;
        end 
    end
end

% Instantiate and Return the Graph Object
G_sys = Graph(sys_verts, sys_edges);

function ConnectX = formatConnectX(ConnectX, class, prop)
    if all(cellfun(@(x) isa(x, class), ConnectX),'all') % Connect X Already given as lists of equivalent GraphVertices or GraphEdges with dominant element in front of list
        return
    elseif all(cellfun(@(x) isa(x, 'double'), ConnectX),'all') % ConnectX specified with component indices in 1st row, Property indices in second row
        num_connections = size(ConnectX,2);
        if size(ConnectX,1) == 3 % If third row specifying dominant component is given
            ConnectX_temp = cell(2, num_connections);
            for c = 1:num_connections
                comps = ConnectX{1,c};
                props = ConnectX{2,c};
                dom_comp = ConnectX{3,c};
                
                i = dom_comp == comps ; % Identify dominant component in first row
                
                new_comps = [comps(i) comps(~i)];
                new_props = [props(i) props(~i)];
                
                ConnectX_temp{1,c} = new_comps;
                ConnectX_temp{2,c} = new_props;
            end
            ConnectX = ConnectX_temp;
        end
        
        vec_lengths = cellfun(@numel, ConnectX);
        assert(all(vec_lengths(1,:)==vec_lengths(2,:)), 'Component and Property Vectors must be the same length for each connection');

        ConnectX_temp = cell(num_connections, 1);
        for c = 1:num_connections
            graphs_i = ConnectX{1,c};
            elems_i = ConnectX{2,c};
            elems = arrayfun(@(c,e) G(c).(prop)(e), graphs_i, elems_i, 'UniformOutput', false); % Arrayfun can't handle heterogenous object arrays so nonuniform output required
            elems = [elems{:}]; % Convert cell array to heterogenous object array
            ConnectX_temp{c} = elems;
        end
        ConnectX = ConnectX_temp;
    else
        error('Invalid ConnectX cell array');
    end
end
end