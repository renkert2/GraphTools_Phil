function [ConnectV] = ExtractVxConn(Graph,ConnectE)
% This function will extract the vertex interconnections for the component
% graph models based on the edge interconnections provided in ConnectE.
% assume no pure vertex connections

% extract all component graphs.
Comp      = vertcat(Graph.Comp);
CompGraph = vertcat(Comp.graph);

% initialized matrix to store vertex interconnection information
% first row:  connected component index
% second row: connected vertex index
% third row:  unused
% fourth row: index of unique interconnection

% Column Pairs: even and odd columns come in pairs. ex: columns 3 and 4
% store info about a single vertex interconnection
% Column Quads: Quad multiple of columns store info about the 2
% vertex interceonnections associated with 1 edge interconnection ex
% 5,6,7, and 8 are two vertex interconnections
V = zeros(4,2*2*size(ConnectE,2));

for i = 1:size(ConnectE,2) % fill first two rows V
    xh1 = CompGraph(ConnectE{1,i}(1)).E(ConnectE{2,i}(1),2); % connected head vertex of component 1
    xt1 = CompGraph(ConnectE{1,i}(1)).E(ConnectE{2,i}(1),1); % connected tail vertex of component 1
    xh2 = CompGraph(ConnectE{1,i}(2)).E(ConnectE{2,i}(2),2); % connected head vertex of component 2
    xt2 = CompGraph(ConnectE{1,i}(2)).E(ConnectE{2,i}(2),1); % connected tail vertex of component 2
        
    V(1,4*i-3:4*i)   = repmat(ConnectE{1,i},1,2); % Graph index
    V(2,4*i-3:4*i-2) = [xh1 xh2]; % head index interconnection
    V(2,4*i-1:4*i-0) = [xt1 xt2]; % tail index interconnection  
end

% loop through V to find instances with the same vertex is shared amongst
% multiple edge interconnections
ind = 1; % indicates the unique interconnection index
for i = 1:2*size(ConnectE,2) % loop through each vertex equivalence
    if V(4,2*i) == 0 % V(4,i) is nonzero if it has been uniquely identified
        V(4,2*i-1:2*i) = ind; % assign a unique identifier to the connection
        idx1 = (V(1, 2*i+1:end) == V(1,2*i-1) & V(2, 2*i+1:end) == V(2,2*i-1)); % find if vertex 1 of a connection is repeated elsewhere
        idx2 = (V(1, 2*i+1:end) == V(1,2*i-0) & V(2, 2*i+1:end) == V(2,2*i-0)); % find if vertex 2 of a connection is repeated elsewhere
        idx = idx1 | idx2; % index of the repeated vertex
        NumCon = reshape(repmat(i+1:2*size(ConnectE,2),2,1),1,2*(2*size(ConnectE,2)-i));
        rep = NumCon.*idx;
        rep(rep == 0) = [];
        if any(idx)
            V(4,[2*rep 2*rep-1]) = ind; % assign the repeated connection index to the repeated connection
        end
        ind = ind + 1; % increase the index number
    end
end

% remove duplicate information in V
V_reduced = flipud(unique(flipud(V)','row')'); %% in improved code, this information should be sufficients for the system generation code

ConnectV = cell(2,ind-1);
% reformat ConnectV as a set
for i = 1:ind-1
    idx = V_reduced(4,:) == i;
    ConnectV(:,i)  = {V_reduced(1,idx); V_reduced(2,idx)};    
end


    

end

