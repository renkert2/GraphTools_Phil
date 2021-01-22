clear all; close all; clc

%% Extract Simulink Models
 
[Comps,ConnectP,G] = autoGraphDefine('ToySystem',1);

%% Plot Component Graphs
figure
for i = 1:length(Comps)
    subplot(ceil(length(Comps)/2),2,i)
    plot(Comps(i).graph,'NodeColor','b','EdgeColor','b');
    title(Comps(i).Name)
end

%% Build System Graph and Model
SystemGraph = Combine(Comps,ConnectP);
SysModel = GraphModel(SystemGraph);

%% Plot System Graph and Model
figure
plot(SysModel,'NodeColor','r','EdgeColor','b','DetailedLabels','States');
