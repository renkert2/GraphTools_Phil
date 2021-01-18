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
[g] = Combine(Graphs, ConnectE, ConnectV, 'CopyEdges', true);
toc
%%
figure
plot(g,'NodeColor','r','EdgeColor','b')



