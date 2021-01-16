clear all; close all; clc

%%
 
[Graph,ConnectE,~] = autoGraphDefine('ToySystem',1);

%%
Comps = [Graph.Comp]';
Graphs = [Comps.graph]';

% plot component graphs
figure
for i = 1:length(Comps)
    subplot(ceil(length(Comps)/2),2,i)
    plot(Comps(i).graph)
    title(Comps(i).Name)
end

Sys = GraphModel(Combine(Graphs,ConnectE));
figure
plot(Sys.graph,'NodeColor','r','EdgeColor','b')
