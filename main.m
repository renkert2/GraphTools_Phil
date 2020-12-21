clear all
close all
clc


%% Load Toy System Connection Set
load ToySystem
figure; plot(G)

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


% plot component graphs
figure
for i = 1:length(Comps)
    subplot(3,2,i)
    plot(Comps(i).graph)
    title(Comps(i).Name)
end


%%
Sys = GenSysGraph(Graphs,ConnectV,ConnectE);
figure
plot(Sys.graph,'NodeColor','r','EdgeColor','b')

% test 2

%% Make a larger graph
% N = 10;
% MLupper = repmat({Sys.M(1:Sys.Nv,:)}, 1, N);
% MLlower = repmat({Sys.M(Sys.Nv+1:end,:)}, 1, N);
% VL = repmat(Sys.Vertices,N,1);
% EL = repmat(Sys.Edges,N,1);
% 
% VL = [VL(arrayfun(@(x) isa(x,'GraphVertex_Internal'),VL));VL(arrayfun(@(x) isa(x,'GraphVertex_External'),VL))];
% EL = [EL(arrayfun(@(x) isa(x,'GraphEdge_Internal'),EL));EL(arrayfun(@(x) isa(x,'GraphEdge_External'),EL))];
% 
% GraphL = Graph([blkdiag(MLupper{:});blkdiag(MLlower{:})],VL,EL);
% SysL = GraphModel(GraphL);
% figure; plot(SysL.graph,'NodeColor','r','EdgeColor','b')
% 
%%
% SymbolicSolver(Sys);


%% 
% [Sys.graph.Vertices(:).Description]'
% tank1.graph.Vertices(2).Description
% tank1.graph.Vertices(2).Description = 'Tank 1 Mass';
% tank1.graph.Vertices(2).Description
% [Sys.graph.Vertices(:).Description]'





% notes:
% WE should remove the external edge class and add it as a vertex property.
% we don't define graphs with external edges part of the edge set so it
% shouldn't be included here.

% added in setaccess = private functionality such that only parameters
% that don't affect graph structure can be edited post system generation


%% 


