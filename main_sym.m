clear all
close all
clc


%% Load Toy System Connection Set
load ToySystem
figure; plot(G)

%% Generate and plot component models
load1 = HeatLoad('Name','Load 1','T_init',symParam('T_init', 16));
HX1   = HeatExchanger('Name','HX 1','T1_init',symParam('T_init', 16));
SJ1   = SplitJunction('Name','J2S1','n_in',2,'n_out',1);
SJ2   = SplitJunction('Name','J1S2','n_in',1,'n_out',2);
tank1 = Tank('Name','Tank 1','T_init',50);
tank2 = Tank('Name','Tank 2','T_init',50);


% don't change the order of the components
Comps = [load1; HX1; SJ1; SJ2; tank1; tank2];
Graphs = [Comps(:).Graph];

% define port connections
ConnectP{1} = [Comps(1).Ports(1) Comps(5).Ports(3)];
ConnectP{2} = [Comps(1).Ports(3) Comps(3).Ports(1)];
ConnectP{3} = [Comps(3).Ports(2) Comps(6).Ports(3)];
ConnectP{4} = [Comps(3).Ports(3) Comps(2).Ports(2)];
ConnectP{5} = [Comps(2).Ports(4) Comps(4).Ports(1)];
ConnectP{6} = [Comps(4).Ports(2) Comps(5).Ports(1)];
ConnectP{7} = [Comps(4).Ports(3) Comps(6).Ports(1)];

% plot component graphs
figure
for i = 1:length(Comps)
    subplot(3,2,i)
    plot(Comps(i).Graph);
    title(Comps(i).Name)
end


%%
SysGraph = Combine(Comps,ConnectP');
Sys = GraphModel(SysGraph);
figure
plot(Sys,'NodeColor','r','EdgeColor','b','DetailedLabels','All');

