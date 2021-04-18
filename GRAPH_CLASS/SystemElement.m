classdef SystemElement < matlab.mixin.Heterogeneous & handle
    %SYSTEMELEMENT Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Name (1,1) string = "Component" % Block Name
        Graph (1,1) Graph
        Ports (:,1) ComponentPort = ComponentPort.empty()
        
        extrinsicProps (:,1) extrinsicProp
        Params (:,1) compParam
    end
    
    methods
        function set.Name(obj, name)
            obj.Name = string(name);
        end
    end
    
    methods (Sealed) 
        function [gSys, extProps] = Combine(C, ConnectP, varargin)
            % Component.Combine Connects components in Component array C according to Port connections defined in ConnectP
            % - ConnectP: Each element of the cell array is a (1xn) Port array  of ports to be connected
            % - varargin optional, passed to Graph.Combine.
            % Returns:
            % - gSys: System Graph
            % - extProps: Combined system extrinsic properties related to the connections
            
            arguments
                C (:,1) SystemElement % Array of elements to be connected in a system
                ConnectP (:,1) cell % Cell array of Component Ports to be connected.
            end
            arguments (Repeating)
                varargin
            end
            
            num_c = size(ConnectP,1);
            ConnectE = {};
            ConnectV = {};
            
            for c = 1:num_c % For each port connection
                ports = ConnectP{c};
                type = ports(1).Type;
                
                assert(isa(ports, 'ComponentPort'), 'Entry %d in ConnectP must be of ComponentPort type',c) 
                assert(all(type == [ports(2:end).Type]), 'Incompatible port types in connection %d', c);
                assert(isCompatible([ports.Domain]), 'Incompatible port domains in connection %d', c);
                
                if type == "EdgeConnection"
                    assert(numel(ports) == 2, 'Edge Connection can only contain two equivalent edges');
                    ConnectE{end+1,1} = [ports.Element];
                elseif type == "VertexConnection"
                    ConnectV{end+1,1} = [ports.Element];
                end 
            end
            
            % Generate System Graph with Combine(G, ConnectE)
            G = [C.Graph];
            gSys = Combine(G, ConnectE, ConnectV, varargin{:});
            
            % Use extrinsicProp.Combine if extProps output argument called
            if nargout == 2
                props = vertcat(C.extrinsicProps);
                extProps = Combine(props);
            end
        end
        
        function setParamVal(obj, param, val)
            obj.Params.setVal(param,val);
        end
    end
end

