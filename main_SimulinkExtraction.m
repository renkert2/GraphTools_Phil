clear all; close all; clc

%%
 
[Graph,ConnectE,~] = autoGraphDefine('ToySystem',1);

%%
Graphs = [vertcat(Graph.Comp).graph];

% plot component graphs
figure
for i = 1:length(Graphs)
    subplot(ceil(length(Graphs)/2),2,i)
    plot(Graphs(i))
    title(Graphs(i).Parent.Name)
end

Sys = GraphModel(Combine(Graphs,ConnectE));
figure
plot(Sys,'NodeColor','r','EdgeColor','b')
