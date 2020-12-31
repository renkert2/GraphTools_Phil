function [Sys, PropMaps] = GenSysGraph_Phil(G, ConnectV, ConnectE) % Input is vector of GraphClass elements
% G - GraphClass Class Array
% ConnectX - Cell vector with first list containing component
% indices; second list containing equivalent properties (edges
% or vertices) corresponding to the components.

arguments
    G (:,1) Graph
    ConnectV cell
    ConnectE cell
end

% Format ConnectV and ConnectE as {:,1} cell array containing list of
% equivalent vertices and edges, respectively
ConnectV = formatConnectX(ConnectV, 'GraphVertex', 'Vertices');
ConnectE = formatConnectX(ConnectE, 'GraphEdge', 'Edges');

% Assign Primary Vertices
primary_vertices = GraphVertex.empty();
for cv = 1:length(ConnectV)
    verts = ConnectV{cv};
    int_vert_flags = arrayfun(@(x) isa(x,'GraphVertex_Internal'));
    n_int_verts = sum(int_vert_flags);
    if n_int_verts == 0
         primary_vertex = verts(1);
    elseif n_int_verts == 1
        primary_vertex = verts(int_vert_flags);
    else
        error('Only one Internal Vertex can be involved in a connection')
    end
    primary_vertices(cv,1) = primary_vertex; 
end


        

% Edge Property Map
[E_map,Ne] = EdgePropMap();

Vertices(Nv) = GraphVertex;
Edges(Ne) = GraphEdge;

for g_i = length(G):-1:1 % Loop backwards so earlier components overide later ones
    Vertices(V_map{g_i}) = G(g_i).Vertices;
    
    Edges_temp = G(g_i).Edges;
    Edges_temp.HeadVertex = 
    Edges_temp.TailVertex = 
    Edges(E_map{g_i}) = G(g_i).Edges;
end

% Call the GraphModel init function to calculate F_inner,etc
Sys.init();

function ConnectX = formatConnectX(ConnectX, class, prop)
    if all(cellfun(@(x) isa(x, class), ConnectX),'all')
        return
    elseif all(cellfun(@(x) isa(x, 'double'), ConnectX),'all')
        vec_lengths = cellfun(@numel, ConnectX);
        assert(all(vec_lengths(1,:)==vec_lengths(2,:)), 'Component and Property Vectors must be the same lenght for each connection');

        num_connections = size(ConnectX,2);
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