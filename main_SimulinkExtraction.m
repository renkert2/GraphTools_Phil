clear all; close all; clc

%%
 
[Comps,ConnectP,G] = autoGraphDefine('ToySystem',1);

%%
Graphs = [Comps.graph];

% plot component graphs
figure
for i = 1:length(Graphs)
    subplot(ceil(length(Graphs)/2),2,i)
    plot(Graphs(i))
    title(Graphs(i).Parent.Name)
end

SystemGraph = Combine(Comps,ConnectP);
SysModel = GraphModel(SystemGraph);
figure
plot(SysModel,'NodeColor','r','EdgeColor','b')
