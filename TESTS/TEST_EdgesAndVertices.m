G = Graph;

V_int(1,1:3) = GraphVertex_Internal;
V_ext(1,1:2) = GraphVertex_External;

V = [V_int V_ext];

G.Vertices = V_int;
G.Vertices = V_ext;
G.Vertices = V;

filter_int = arrayfun(@(x) isa(x,'GraphVertex_Internal'),V); % Find which vertices in V are external
filter_ext = arrayfun(@(x) isa(x,'GraphVertex_External'),V); % Find which vertices in V are external

G.Vertices.Description % Get the description of all vertices;
G.Vertices(filter_int).Capacitance % Get Capacitance of all Internal vertices
% G.Vertices.Capacitance % throws an error because G.Vertices only has the
% 'Description' property in common.
