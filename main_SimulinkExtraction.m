clear all; close all; clc

%%
[Graph,~] = autoGraphDefine('ToySystem',1);

%%
ConnectE = ExtractExConn(Graph);
ConnectV = ExtractVxConn(Graph,ConnectE);

Comps = [Graph.Comp]';
Graphs = [Comps.graph]';

% plot component graphs
figure
for i = 1:length(Comps)
    subplot(3,2,i)
    plot(Comps(i).graph)
    title(Comps(i).Name)
end


Sys = GenSysGraph(Graphs,ConnectV,ConnectE);
figure
plot(Sys.graph,'NodeColor','r','EdgeColor','b')