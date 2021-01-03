clear all
close all
clc


%% Load Toy System Connection Set
load ToySystem
% figure; plot(G)

%% Generate and plot component models

HX1   = HeatExchanger('Name','HX 1','T1_init',16);
load1 = HeatLoad('Name','Load 1','T_init',16);
SJ1   = SplitJunction('Name','J2S1','n_in',2,'n_out',1);
SJ2   = SplitJunction('Name','J1S2','n_in',1,'n_out',2);
tank1 = Tank('Name','Tank 1','T_init',50);
tank2 = Tank('Name','Tank 2','T_init',50);


% don't change the order of the components
Comps = [HX1; load1; SJ1; SJ2; tank1; tank2];
Graphs = [Comps(:).graph];

%%
tic
[Sys_graph_phil, ConnectE_mod, ConnectV_mod] = GenSysGraph_Phil(Graphs, ConnectE, ConnectV, 'CopyEdges', true);
new_time = toc

tic 
[Sys_graph_phil] = GenSysGraph_Phil(Graphs, ConnectE_mod, ConnectV_mod(1:6), 'CopyEdges', true);
new_time_2 = toc

tic
[Sys_graph_phil] = Combine(Graphs, ConnectE_mod, ConnectV_mod(1:6), 'CopyEdges', true);
new_time_3 = toc

tic
Sys_graph_old = GenSysGraph(Graphs, ConnectV, ConnectE);
old_time = toc

%% Ensure both graphs are equivalent
verts = [vertcat(Sys_graph_phil.Vertices) vertcat(Sys_graph_old.Vertices)];
edges = [vertcat(Sys_graph_phil.Edges) vertcat(Sys_graph_old.Edges)];

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


%%
figure
plot(Sys_graph_phil,'NodeColor','r','EdgeColor','b')

figure
plot(Sys_graph_old,'NodeColor','r','EdgeColor','b')

Sys = GraphModel(Sys_graph_old);



