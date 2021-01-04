function graph_connections()
%GRAPH_CONNECTIONS initializes Graph ID block
%
%   GRAPH_CONNECTIONS runs in the mask initalization of each Graph ID 
%   block. The function updates the mask parameters of the initializing 
%   mask. The function executes fully when SimulationStatus is updating,
%   otherwise the function exits early.
%
%   Written by: Christopher T. Aksland
%   

%% CHECK MODEL STATUS -- only run the code if model is updating
if ~strcmp(get_param(bdroot, 'SimulationStatus'),'updating')
    return
end

%% INITIALIZATION
downType = '';  downPath = '';
upType = '';    upPath = '';

%% GRAB MODEL HANDLES AND PARAMETERS
H               = get_param(gcb,'handle');  % Get block handle
Hblock_current  = get(H);                   % Parse out handles structure
H_parent        = get_param(Hblock_current.Parent,'handle');
Hblock_parent   = get(H_parent);

% Mask values for the currently updating mask
currMaskValues = get_param(H,'MaskValues');

% loop through ports of the parent block
for i=1:numel(Hblock_parent.PortConnectivity)
    
    % TRY fails if upstream or downstream block has no children named 
    % Graph ID -- no need for conditional statement
    try
        if ~isempty(Hblock_parent.PortConnectivity(i).DstBlock)
            % downstream block parameters
            hblock_dwn = get(Hblock_parent.PortConnectivity(i).DstBlock);
            
            % find blocks named Graph ID  in downstream block - get param
            down = get(get_param([hblock_dwn.Path,'/',hblock_dwn.Name,'/Graph ID'],'handle'));

            % concatenate the downstream block type & path
            downType = [downType,down.MaskValues{2},','];
            downPath = [downPath,down.Path,','];

        elseif ~isempty(Hblock_parent.PortConnectivity(i).SrcBlock)
            % upstream block parameters
            hblock_up = get(Hblock_parent.PortConnectivity(i).SrcBlock);

            % find blocks named OMEGA Graph in upstream block - get param
            up = get(get_param([hblock_up.Path,'/',hblock_up.Name,'/Graph ID'],'handle'));

            % concatenate the upstream block type & path
            upType = [upType,up.MaskValues{2},','];
            upPath = [upPath,up.Path,','];
            
        end
    catch % no conditions for the CATCH
    end
    
end

% update current mask values with upstream/downstream type and path
currMaskValues{3} = downType(1:end-1);
currMaskValues{4} = downPath(1:end-1);
currMaskValues{5} = upType(1:end-1);
currMaskValues{6} = upPath(1:end-1);
set_param(H,'MaskValues',currMaskValues);