clear all
close all
clc


%% Load Toy System Connection Set
load ToySystem

%% Generate and plot component models
HX1   = HeatExchanger('Name','HX 1','T1_init',16);
load1 = HeatLoad('Name','Load 1','T_init',16);
SJ1   = SplitJunction('Name','J2S1','n_in',2,'n_out',1);
SJ2   = SplitJunction('Name','J1S2','n_in',1,'n_out',2);
tank1 = Tank('Name','Tank 1','T_init',50);
tank2 = Tank('Name','Tank 2','T_init',50);

% don't change the order of the components
Comps = [HX1; load1; SJ1; SJ2; tank1; tank2];
Graphs = [Comps.graph];

%% Add Ports to Components - Pulled from ConnectE
% Component 1: HX 1
p1(1) = ComponentPort('Element', Graphs(1).Edges(2));
p1(2) = ComponentPort('Element', Graphs(1).Edges(1));

p1(3) = ComponentPort('Element', Graphs(1).Vertices(5));
p1(4) = ComponentPort('Element', Graphs(1).Vertices(1));
p1(5) = ComponentPort('Element', Graphs(1).Vertices(3));
% Component 2: Load 1
p2(1) = ComponentPort('Element', Graphs(2).Edges(2));
p2(2) = ComponentPort('Element', Graphs(2).Edges(1));

p2(3) = ComponentPort('Element', Graphs(2).Vertices(3));
p2(4) = ComponentPort('Element', Graphs(2).Vertices(1));
p2(5) = ComponentPort('Element', Graphs(2).Vertices(2));
% Component 3: SJ 1
p3(1) = ComponentPort('Element', Graphs(3).Edges(1));
p3(2) = ComponentPort('Element', Graphs(3).Edges(3));
p3(3) = ComponentPort('Element', Graphs(3).Edges(2));

p3(4) = ComponentPort('Element', Graphs(3).Vertices(4));
p3(5) = ComponentPort('Element', Graphs(3).Vertices(1));
p3(6) = ComponentPort('Element', Graphs(3).Vertices(2));
p3(7) = ComponentPort('Element', Graphs(3).Vertices(3));
% Component 4: SJ 2
p4(1) = ComponentPort('Element', Graphs(4).Edges(1));
p4(2) = ComponentPort('Element', Graphs(4).Edges(2));
p4(3) = ComponentPort('Element', Graphs(4).Edges(3));

p4(4) = ComponentPort('Element', Graphs(4).Vertices(1));
p4(5) = ComponentPort('Element', Graphs(4).Vertices(2));
p4(6) = ComponentPort('Element', Graphs(4).Vertices(3));
p4(7) = ComponentPort('Element', Graphs(4).Vertices(4));
% Component 5: Tank 1
p5(1) = ComponentPort('Element', Graphs(5).Edges(1));
p5(2) = ComponentPort('Element', Graphs(5).Edges(2));

p5(3) = ComponentPort('Element', Graphs(5).Vertices(3));
p5(4) = ComponentPort('Element', Graphs(5).Vertices(4));
p5(5) = ComponentPort('Element', Graphs(5).Vertices(1));
% Component 6: Tank 2
p6(1) = ComponentPort('Element', Graphs(6).Edges(1));
p6(2) = ComponentPort('Element', Graphs(6).Edges(2));

p6(3) = ComponentPort('Element', Graphs(6).Vertices(3));
p6(4) = ComponentPort('Element', Graphs(6).Vertices(4));
p6(5) = ComponentPort('Element', Graphs(6).Vertices(1));

% Assign ports to Components
Comps(1).Ports = p1;
Comps(2).Ports = p2;
Comps(3).Ports = p3;
Comps(4).Ports = p4;
Comps(5).Ports = p5;
Comps(6).Ports = p6;

%% Create ConnectP
ConnectP = {[p1(1) p4(1)];
        [p2(1) p3(1)];
        [p3(2) p1(2)];
        [p4(2) p6(1)];
        [p4(3) p5(1)];
        [p5(2) p2(2)];
        [p6(2) p3(3)];
        [p4(4) p5(3) p6(3) p1(3)];
        [p1(4) p4(5) p3(4)];
        [p3(5) p1(5) p2(3) p6(4)];
        [p2(4) p3(6) p5(4)];
        [p6(5) p3(7) p4(6)];
        [p5(5) p2(5) p4(7)]};
%% Generate Gsys from Graph.Combine
Gsys = Combine(Comps, ConnectP, 'CopyEdges', true);
Gsys_old = Combine(Graphs, ConnectE, ConnectV, 'CopyEdges', true);

%% Ensure both graphs are equivalent
verts = [vertcat(Gsys.Vertices) vertcat(Gsys_old.Vertices)];
edges = [vertcat(Gsys.Edges) vertcat(Gsys_old.Edges)];

nv = size(verts,1);
vert_grid = zeros(nv);

for i = 1:nv
    for j = 1:nv
        vert_grid(i,j) = isequal(verts(i,1), verts(j, 2));
    end
end

ne = size(edges,1);
edge_grid = zeros(ne);

for i = 1:ne
    for j = 1:ne
        edge_grid(i,j) = isequal(edges(i,1), edges(j, 2));
    end
end

assert(sum(vert_grid,'all') == nv, 'Vertices do not match!')
assert(sum(edge_grid,'all') == ne, 'Edges do not match!')