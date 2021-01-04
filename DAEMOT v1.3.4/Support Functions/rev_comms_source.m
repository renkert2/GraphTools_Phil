function rev_comms_source(n_flows)

h      = get_param(gcb,'handle');   % Get block handle
hblock_current = get(h);     % Parse out handles structure
h_parent = get_param(hblock_current.Parent,'handle');
hblock_parent = get(h_parent);

%%% downstream block goto tag
n_inPorts = hblock_parent.Ports(1);

if isempty(hblock_parent.PortConnectivity(n_inPorts+n_flows,1).DstBlock)
    hblock_goto.GotoTag = 'NULL';
else
    for i=1:size(hblock_parent.PortConnectivity(n_inPorts+n_flows,1).DstBlock,2)
        hblock_dwn = get(hblock_parent.PortConnectivity(n_inPorts+n_flows,1).DstBlock(i));
        prt = hblock_parent.PortConnectivity(n_inPorts+n_flows,1).DstPort(i);
        if strcmp(hblock_dwn.BlockType,'Goto')
            % get a list of the From blocks
            froms = find_system(bdroot, 'BlockType', 'From');
            % get a list of the From tags
            fromTags = get_param(froms, 'GotoTag');
            gotoTag = hblock_dwn.GotoTag;
            % get the indices of gotoTags that do not have fromTags
            indices = find(strcmp(fromTags,gotoTag));
            % get the first corresponding from tag (should only be one)
            from = froms(indices(1));
            from =  cell2mat(from);
            h_from = get_param(from,'handle');
            hblock_parent = get(h_from);
            hblock_dwn = get(hblock_parent.PortConnectivity(1,1).DstBlock(1));
            for j = 1:max(size(hblock_dwn.PortConnectivity))
                if hblock_dwn.PortConnectivity(j).SrcBlock == h_from
                    prt = hblock_dwn.PortConnectivity(j).Type;
                end
            end
            prt = str2num(prt) - 1;
        end
        %     str = [bdroot,'/',hblock_dwn.Name,'/Sink',num2str(prt+1),'/Goto'];
        str = [hblock_dwn.Path,'/',hblock_dwn.Name,'/Sink',num2str(prt+1),'/Goto'];
        try
            h = get_param(str,'handle');
            hblock_goto = get(h);
        catch
            %error('test')
        end
    end
end

%%% local goto block
str = [gcb,'/From'];
set_param(str,'GotoTag',hblock_goto.GotoTag);