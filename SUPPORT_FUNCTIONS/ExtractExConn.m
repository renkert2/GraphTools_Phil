function [ConnectE] = ExtractExConn(Graph)
% This function will extract the edge interconnections for the component
% graph models based on the structure of a Simulink Model

% ConnectE(:,i)  = {[x y]    ;[E_x_m E_y_n]; Z};
% for x,y,... in {1,2,...,N} for N graphs
% E indicates Edge
% n,m indicate n-th and m-th edges of graphs x and y respectively
% Z is an element of {x,y}
% the first cell indicates which graphs are interconnected
% the second cell indicates which edges of the graphs are equivalent
% the third cell indicates which edge type is resultant in the equivalency
%   (for now, the third cell is just the first graph in the connection)

%%%%%%%% Here, "tail and head vertex" refer to the graph representation of
%%%%%%%% the block model, not the graph model

% first, remove the elements of graph that do not have associated component
% models (ie. remove sink and source elements).
Comp = vertcat(Graph.Comp);
nComp = length(Comp); % number of elements with components
Graph(nComp+1:end) = []; % remove elements without components
for i = 1:nComp
    Graph(i).DownVertex(Graph(i).DownVertex > nComp) = 0; % replace up and downstream vertex idicies of sinks/sources to 0
    Graph(i).UpVertex(Graph(i).UpVertex > nComp) = 0;  % replace up and downstream vertex idicies of sinks/sources to 0
    Graph(i).Port = [Graph(i).UpVertex, Graph(i).DownVertex]; % ports is all upstream and downsteam vertices

end

% initialize variables
E = zeros(2); % Edge matrix for the block diagram
ind = 1; 
ConnectE = cell(3,size(E,1));

CompGraph = vertcat(Comp.graph); % list of all component graphs
for i = 1:numel(Graph) % loop through Component Graphs to develop edge connection sets 
    for j = 1:numel(Graph(i).DownVertex) % loop through all downstream connections
        xT = i; % tail vertex of block digram
        xH = Graph(i).DownVertex(j); % head vertex of block diagram
        
        if xH ~= 0 % incase there is an unattached vertex ie. if the connection involves a source or sink, skip it
            % idxH and idxT are Port indices
            % ex: this connection connects port # 'idxH' of the head block 
            % model and port # 'idxT' of the tail block block
            idxH = find(xT == Graph(xH).Port); % the index of the head vertex port connected to the tail vertex. 
            idxT = find(xH == Graph(xT).Port); % the index of the tail vertex port connected to the head vertex
            
            portT = [CompGraph(xT).InternalEdges(:).Port]; % list of ports associated with edges for tail vertex
            portH = [CompGraph(xH).InternalEdges(:).Port]; % list of ports associated with edges for head vertex
            
            % figure out which graph model edges are associated with each block model port 
            eConn = [find(idxT == portT) find(idxH == portH)]; %[tail vertex port associated edge , head vertex port associated edge]
            
            E(ind,1) = xT; % edge matrix tail
            E(ind,2) = xH; % edge matrix head
            % build ConnectE. Note the 
            ConnectE(:,ind) = {E(ind,:); eConn ;E(ind,1)};
            ind = ind + 1;
        end
    end
end

% figure
% G = digraph(E(:,1),E(:,2));
% h = plot(G);
% labelnode(h,1:numel(Graph),{Graph(:).Name})



end

