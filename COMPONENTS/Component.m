classdef Component < matlab.mixin.Heterogeneous & handle
    %COMPONENT Super class to be inherited by Components (i.e.
    %Battery, Motor, Heat Exchanger, etc...
    %   Detailed explanation goes here
    
    properties %(SetAccess = protected)        
        % Block Name
        Name string = string.empty()
        graph Graph = Graph.empty()
        Ports ComponentPort = ComponentPort.empty()
    end
    
    methods
        function obj = Component(varargin)
            if nargin == 1
                if isstruct(varargin{1})
                    for i = 1:numel(varargin{1})
                        try
                            obj.(varargin{1}(i).Name) = varargin{1}(i).Value;
                        catch
                            % eventually update this to indicate that no
                            % property for this class exists.
                        end
                    end
                else
                    error('Components must be defined using a structure of Name/Value fields or Name-Value pairs')
                end
            elseif nargin > 1
                my_inputparser(obj,varargin{:}); % input parser component models
            end
            obj.init(); % I don't know why we need this and can't just call ConstructGraph - CTA

        end
        
        function set.Name(obj, name)
            obj.Name = string(name);
        end
        
        function init(obj)
            obj.ConstructGraph();
            obj.DefineChildren();
        end
       
    end
    
    methods (Sealed)
         
        function gSys = Combine(C, ConnectP, varargin)
            arguments
                C (:,1) Component % Array of components to be connected in a system
                ConnectP (:,1) cell % vector of Component Ports to be connected.  Connections along dimension 1, equivalent ports along dimension 2
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
                domain = ports(1).Domain;
                
                assert(isa(ports, 'ComponentPort'), 'Entry %d in ConnectP must be of ComponentPort type',c) 
                assert(all(type == [ports(2:end).Type]), 'Incompatible port types in connection %d', c);
                assert(all(domain == [ports(2:end).Domain]), 'Incompatible port types in connection %d', c);
                
                if type == 1 % Type 1 Connection
                    assert(numel(ports) == 2, 'Type 1 Connection can only contain two edges');
                    ConnectE{end+1,1} = [ports.Element];
                elseif type == 2 % Type 2 Connection
                    ConnectV{end+1,1} = [ports.Element];
                end
                
                
            end
            
            % Construct G
            G = [C.graph];
            
            % Generate System Graph with Combine(G, ConnectE)
            gSys = Combine(G, ConnectE, ConnectV, varargin{:});
            
        end
    end
    
    methods (Access = protected)
        function ConstructGraph(obj)
            g = DefineGraph(obj);
            obj.graph = g;
        end
        
        function DefineChildren(obj)
            try
                obj.graph.Parent = obj;
                graph_children = ["Vertices", "Edges", "Inputs"];
                for child = graph_children
                    for i = 1:numel(obj.graph.(child))
                        obj.graph.(child)(i).Parent = obj;
                    end
                end
            catch
                warning('Error defining component as parent object')
            end
        end
    end
    
    methods (Abstract, Access = protected)
        DefineGraph(p)       
    end
end

