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

% plot component graphs
figure
for i = 1:length(Comps)
    subplot(3,2,i)
    plot(Comps(i).Model)
    title(Comps(i).Name)
end

Graphs = [Comps(:).Model];

%%
Sys = GenSysGraph(Graphs,ConnectV,ConnectE);
figure
plot(Sys,'NodeColor','r','EdgeColor','b')
Sys.MakeMatrices()
% test 2


%% 
% [Sys.Graph.Vertices(:).Description]'
% tank1.Model.Graph.Vertices(2).Description
% tank1.Model.Graph.Vertices(2).Description = 'Tank 1 Mass';
% tank1.Model.Graph.Vertices(2).Description
% [Sys.Graph.Vertices(:).Description]'





% notes:
% WE should remove the external edge class and add it as a vertex property.
% we don't define graphs with external edges part of the edge set so it
% shouldn't be included here.

% added in setaccess = private functionality such that only parameters
% that don't affect graph structure can be edited post system generation


%% 


