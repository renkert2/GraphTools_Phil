function [COMPS,ConnectP,G] = autoGraphDefine(SYS,PLOTGRAPH)
%AUTOGRAPHDEFINE generates directed graph of DAEMOT system
%
%   [G, GRAPH] = AUTOGRAPHDEFINE(SYS) parses the Simulink model SYS for all 
%   Ewok blocks with an Graph ID block. Using these block, 
%   AUTOGRAPHDEFINE generates a directed graph G and the full data 
%   structure GRAPH
%
%   [G, GRAPH] = AUTOGRAPHDEFINE(SYS,PLOTGRAPH) parses the Simulink model 
%   SYS for all DAEMOT blocks with an Graph ID block, generates, and plots 
%   the directed graph G if PLOTGRAPH is greater than 0.
%
%   Written by: Christopher T. Aksland
%  

% step 1: update the Simulink model to force all Graph ID blocks to
% update the upstream/downstream connections
set_param(SYS, 'SimulationCommand', 'update')

% step 2: find all Graph ID blocks within the system
Objects = find_system(SYS,'LookUnderMasks','on','FollowLinks','on','Name','Graph ID');
                      
% step 3: parse each Graph ID block and pull out parameters
GRAPH = struct;     % init Graph as a struct

% loop over all the found Graph ID blocks
for i = 1:numel(Objects)
    hObj = get_param(Objects{i},'handle');      % Graph ID block handle
    objParam = get(hObj);                       % Graph ID block params
    Name = split(objParam.Path,'/');            % Block name
    GRAPH(i).Vertex = i;                        % vertex number
    GRAPH(i).Name = Name{end};                  % vertex name
    GRAPH(i).Type = convertCharsToStrings(objParam.MaskValues{2});     % vertex type
    GRAPH(i).Path = {objParam.Path};              % Block path
    GRAPH(i).downType = objParam.MaskValues{3}; % downstream block type
    GRAPH(i).downPath = split(objParam.MaskValues{4},','); % down blk path
    GRAPH(i).upType = objParam.MaskValues{5};   % upstream block type
    GRAPH(i).upPath = split(objParam.MaskValues{6},',');   % ups blk path
    GRAPH(i).Special = objParam.MaskValues{7};  % special notes
    
    % build the component graph models
    try
        params = get(getSimulinkBlockHandle(GRAPH(i).Path),'MaskWSVariables');
        params(end+1).Name = 'Name'; params(end).Value = GRAPH(i).Name;
        GRAPH(i).Comp = feval(GRAPH(i).Type,params);
    catch
        
    end
        
end

% remove the Special field
GRAPH = rmfield(GRAPH, 'Special');

%% Step 3.5.5 Reorder components
% reorder components such that physical components are ordered before sinks
% and sources 
idx = cellfun(@isempty,{GRAPH(:).Comp}); %index location of physical elements
GRAPH = [GRAPH(~idx), GRAPH(idx)]; % reorder the components
for i = 1:numel(GRAPH)
    GRAPH(i).Vertex = i; % relabel the Vertex property of the GRAPH structure
end


%% step 4: loop through each Object and compare path to the other object
% upstream and downstream paths -- find adjacent vertices
% the path property is the unique identifier of each element
for i = 1:numel(GRAPH)
    % current component, upstream, & downstream paths
    currDownPath = GRAPH(i).downPath;
    currUpPath = GRAPH(i).upPath;
    
    % intialize parameters
    GRAPH(i).DownVertex = [];
    GRAPH(i).UpVertex = [];
    
    for j = 1:numel(GRAPH)  % loop through each element of Graph
        if i ~= j           % skip comparing an element to itself
            for k = 1:numel(currDownPath) % loop over multiple paths
                % if downstream path of current block matches the path of
                % another block, then the components are connected and j is
                % the downstream vertex number 
                if nnz(strcmp(currDownPath{k},vertcat(GRAPH(j).Path{:})))
                    GRAPH(i).DownVertex(k) = j;
                end
            end
            for k = 1:numel(currUpPath)
                % if upstream path of current block matches the path of
                % another block, then the components are connected and j is
                % the upstream vertex number 
                if nnz(strcmp(currUpPath{k},vertcat(GRAPH(j).Path{:})))
                    GRAPH(i).UpVertex(k) = j;
                end
            end
        end
    end
end

%% step 5: define the edge matrix, E
E = zeros(2); ind = 1;  % initialize

% loop over each element in Graph
for i = 1:numel(GRAPH)
    % loop over the number of downstream vertices for the i-th element
    for j = 1:numel(GRAPH(i).DownVertex)
        E(ind,1) = GRAPH(i).Vertex; % edge originates at the i-th element
        E(ind,2) = GRAPH(i).DownVertex(j); % edge ends at downstream vertex
        ind = 1 + ind;
    end
end

% generate the directional graph
G = digraph(E(:,1),E(:,2),ones(numel(E)/2,1),{GRAPH.Name});   


%% PLOT THE GRAPH
if PLOTGRAPH
    figure; p = plot(G,'Layout','force');
%     set(gcf,'Position',[545 252 1228 559]);
    
    names = p.NodeLabel; % same node labels for later
    p.NodeLabel = {}; p.NodeColor = 'k'; p.MarkerSize = 8;  % marker properties
    p.EdgeColor = 'k'; p.LineWidth = 1.5; p.ArrowSize = 15; % edge properties
    
    % Custom labels
    hold on
    set(p, 'DefaultTextInterpreter', 'none')
    for i=1:length(names)
        text(p.XData(i), p.YData(i), ['  ',names(i)], 'Color', 'k', 'FontSize', 10);
    end
    hold off
end

%% Get Component Connection Set
% edge matrix for internal connections
E_int = E((E(:,1)<=sum(~idx) & E(:,2)<=sum(~idx)),:);
% add Port information
Comp = vertcat(GRAPH.Comp);
nComp = length(Comp); % number of elements with components
for i = 1:nComp
    GRAPH(i).Port = [GRAPH(i).UpVertex, GRAPH(i).DownVertex]; % ports is all upstream and downsteam vertices
end

ConnectP = cell(1,size(E_int,1));
for i = 1:length(ConnectP)
    port1 = find(GRAPH(E_int(i,1)).Port==E_int(i,2)); % port of the first component connected to second component
    port2 = find(GRAPH(E_int(i,2)).Port==E_int(i,1)); % port of the second component connected to first component
    ConnectP{i} = [GRAPH(E_int(i,1)).Comp.Ports(port1) GRAPH(E_int(i,2)).Comp.Ports(port2)]; 
end

%% 
COMPS = vertcat(GRAPH.Comp);
% ConnectE = ExtractExConn(GRAPH);

