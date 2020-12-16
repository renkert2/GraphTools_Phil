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
tank2 = Tank('Name','Tank 1','T_init',50);

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
plot(Sys)
% test 2