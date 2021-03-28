classdef Component < matlab.mixin.Heterogeneous & handle
    % Component is a super class for all specifc components (Ex: tank, 
    % battery, etc) in the Graph Modeling Toolbox. The graph property of
    % different components can be connected to generate a system model.
    % Instatiate an empty object, use an input parser, or use a strcuture array
    % with "Name" and "Value" fields, or use a single structure with 
    % fields corresponding to parameters. Valid Name and Value pairs are
    % dictated by the component subclasses.
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Contributors: Christopher T. Aksland and Phil Renkert
    % Association: University of Illionis at Urbana-Champaign
    % Contact: aksland2@illinois.edu and renkert2@illinois.edu
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % potential improvements:
    % Move SymParam and extrinsicProps outside of GraphClass Core
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    properties       
        Name (1,1) string = "Component" % Block Name
        Graph (1,1) Graph
        Ports (:,1) ComponentPort = ComponentPort.empty()
        
        extrinsicProps (:,1) extrinsicProp
        SymParams SymParams {mustBeScalarOrEmpty}
    end
    
    methods
        function obj = Component(varargin)
            if nargin == 1
                if isstruct(varargin{1}) % create a component using a strucutre
                    
                    if isequal(fields(varargin{1}),{'Name';'Value'}) % Structure passed as struct array with property 'Name' and 'Value' fields, e.x. s(1).Name = 'R', s(1).Value = 1
                        for i = 1:numel(varargin{1})
                            try
                                obj.(varargin{1}(i).Name) = varargin{1}(i).Value;
                            catch
                                % eventually update this to indicate that no
                                % property for this class exists.
                            end
                        end
                        
                    else
                        fnames = fieldnames(varargin{1}); % Single structure passed with each field corresponding to a property, e.x. s.R = 1
                        for i = 1:numel(fnames)
                            try
                                obj.(fnames{i}) = varargin{1}.(fnames{i});
                            end
                        end
                    end
                else
                    error('Components must be defined using a struct array with Name/Value fields or a single struct with fields corresponding to properties')
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
            % Dive into obj.Graph and set the Parent property of all objects with the Parent property
            
            try
                obj.Graph.Parent = obj;
                graph_children = ["Vertices", "Edges", "Inputs", "Outputs"];
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
            % Pass symbolic parameters defined in the Component properties 
            % to the Graph
            props = properties(obj);
            sp_cell = {};
            for i = 1:numel(props)
                prop = obj.(props{i});
                if isa(prop, 'symParam')
                    sp_cell{end+1} = prop;
                end
            end
            sym_params = SymParams(sp_cell);
            obj.SymParams = sym_params;
            obj.Graph.SymParams = sym_params;
        end          
    end
    
        
    methods (Access = protected)        
        function init(obj)
            % Placeholder - Optionally modified in subclasses
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
                C (:,1) Component % Array of components to be connected in a system
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
            obj.SymParams.setVal(param,val);
        end
        
        function comp_arrays = Replicate(obj_array, N)
            % Replicate is used to create N independent copies of all the components in obj_array
            % Returns comp_arrays cell array containing the duplicate components.  
            % - comp_arrays{i} is an 1xN Component array of copies corresponding to obj_array(i)
            
            comp_arrays = cell(size(obj_array));
            for i = 1:numel(obj_array)
                comp = obj_array(i);
                unique_props = setdiff(properties(comp), properties('Component')); % Get properties unique to the specific component, i.e. remove properties in the Component superclass
                prop_struct = struct();
                for j = 1:numel(unique_props)
                    prop_struct.(unique_props{j}) = comp.(unique_props{j});
                end
                
                constructFun = str2func(class(comp)); % Get the class Constructor handle specific to the component.
                
                for j = 1:N
                    prop_struct.Name = join([comp.Name,string(j)]);
                    comp_array(j) = constructFun(prop_struct); % Call the constructor N times to get N independent objects with identical property values.
                end
                
                comp_arrays{i} = comp_array;
            end
        end
    end
    
    methods (Abstract, Access = protected)
        DefineComponent(obj) % Required in subclasses
    end
end

