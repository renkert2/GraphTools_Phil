classdef Component < matlab.mixin.Heterogeneous & handle
    %COMPONENT Super class to be inherited by Components (i.e.
    %Battery, Motor, Heat Exchanger, etc...
    %   Detailed explanation goes here
    
    properties %(SetAccess = protected)        
        % Block Name
        Name string = "Component"
        Graph Graph = Graph.empty()
        Ports ComponentPort = ComponentPort.empty()
        
        extrinsicProps (:,1) extrinsicProp
    end
    
    methods
        function obj = Component(varargin)
            if nargin == 1
                if isstruct(varargin{1})
                    if isequal(fields(varargin{1}),{'Name';'Value'})
                        for i = 1:numel(varargin{1})
                            try
                                obj.(varargin{1}(i).Name) = varargin{1}(i).Value;
                            catch
                                % eventually update this to indicate that no
                                % property for this class exists.
                            end
                        end
                    else
                        fnames = fieldnames(varargin{1});
                        for i = 1:numel(fnames)
                            try
                                obj.(fnames{i}) = varargin{1}.(fnames{i});
                            end
                        end
                    end
                else
                    error('Components must be defined using a structure of Name/Value fields or Name-Value pairs')
                end
            elseif nargin > 1
                my_inputparser(obj,varargin{:}); % input parser component models
            end
            obj.init_super();

        end
        
        function set.Name(obj, name)
            obj.Name = string(name);
        end
    end
    
    methods (Sealed)
        function init_super(obj)
            obj.init(); % Concrete method that can be overriden by subclasses
            obj.DefineComponent();
            obj.DefineChildren();
            obj.DefineSymParams();
        end
    end
    
    methods (Sealed, Access = protected)
        function DefineChildren(obj)
            try
                obj.Graph.Parent = obj;
                graph_children = ["Vertices", "Edges", "Inputs","Outputs"];
                for child = graph_children
                    for i = 1:numel(obj.Graph.(child))
                        obj.Graph.(child)(i).Parent = obj;
                    end
                end
            catch
                warning('Error defining component as parent object')
            end
        end
        
        function DefineSymParams(obj)
            props = properties(obj);
            
            sym_params = sym.empty();
            sym_params_vals = [];
            
            for i = 1:numel(props)
                prop = obj.(props{i});
                if isa(prop, 'symParam')
                    sym_params(end+1,1) = prop;
                    sym_params_vals(end+1,1) = double(prop);
                end
            end
            
            obj.Graph.SymParams = sym_params;
            obj.Graph.SymParams_Vals = sym_params_vals;
        end          
    end
    
        
    methods (Access = protected)        
        function init(obj)
            % Placeholder - Optionally modified in subclasses
        end
    end
    
    methods (Sealed) 
        function [gSys, extProps] = Combine(C, ConnectP, varargin)
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
                
                assert(isa(ports, 'ComponentPort'), 'Entry %d in ConnectP must be of ComponentPort type',c) 
                assert(all(type == [ports(2:end).Type]), 'Incompatible port types in connection %d', c);
                assert(isCompatible([ports.Domain]), 'Incompatible port domains in connection %d', c);
                
                if type == 1 % Type 1 Connection
                    assert(numel(ports) == 2, 'Type 1 Connection can only contain two edges');
                    ConnectE{end+1,1} = [ports.Element];
                elseif type == 2 % Type 2 Connection
                    ConnectV{end+1,1} = [ports.Element];
                end 
            end
            
            % Generate System Graph with Combine(G, ConnectE)
            G = [C.Graph];
            gSys = Combine(G, ConnectE, ConnectV, varargin{:});
            
            if nargout == 2
                props = vertcat(C.extrinsicProps);
                extProps = Combine(props);
            end
        end
        
        function comp_arrays = Replicate(obj_array, N)
            comp_arrays = cell(size(obj_array));
            for i = 1:numel(obj_array)
                comp = obj_array(i);
                unique_props = setdiff(properties(comp), properties('Component'));
                prop_struct = struct();
                for j = 1:numel(unique_props)
                    prop_struct.(unique_props{j}) = comp.(unique_props{j});
                end
                
                constructFun = str2func(class(comp)); % Get the class constructor specific to the component
                
                for j = 1:N
                    prop_struct.Name = join([comp.Name,string(j)]);
                    comp_array(j) = constructFun(prop_struct); % Call the constructor N times to get N independent objects.
                end
                
                comp_arrays{i} = comp_array;
            end
        end
    end
    
    methods (Abstract, Access = protected)
        DefineComponent(obj) % Required in subclasses
    end
end

