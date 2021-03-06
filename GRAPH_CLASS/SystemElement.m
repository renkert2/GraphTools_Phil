classdef SystemElement < matlab.mixin.Heterogeneous & matlab.mixin.Copyable & handle
    %SYSTEMELEMENT Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Name (1,1) string = "Component" % Block Name
        Graph {mustBeScalarOrEmpty} = Graph.empty()
        Ports (:,1) ComponentPort = ComponentPort.empty()
        
        Params (:,1) compParam
    end
    
    methods
        function obj = SystemElement(varargin)
            if nargin > 0
                my_inputparser(obj, varargin{:});
            end
        end
    end
    
    methods (Sealed) 
        function [Sys, gSys, params] = Combine(C, ConnectP, Sys, varargin)
            % Component.Combine Connects components in Component array C according to Port connections defined in ConnectP
            % - ConnectP: Each element of the cell array is a (1xn) Port array  of ports to be connected
            % - varargin optional, passed to Graph.Combine.
            % Returns:
            % - gSys: System Graph
            % - extProps: Combined system extrinsic properties related to the connections
            
            arguments
                C (:,1) SystemElement % Array of elements to be connected in a system
                ConnectP (:,1) cell % Cell array of Component Ports to be connected.
                Sys SystemElement = SystemElement.empty()
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
            
            if isempty(Sys)
                Sys = System();
            end
            
            Sys.Graph = gSys;
            Sys.Components = C;
            Sys.ConnectP = ConnectP;
            setCombinedParams(Sys);
        end
        
        function params = setCombinedParams(obj,opts)
            arguments
                obj
                opts.FilterDirectDescendentExProps logical = false 
            end
            
            existing_params = obj.Params;
            C = obj.Components;
            params = vertcat(C.Params);
            if ~isempty(params) % Params have components
                exprop_filter = isaArrayFun(params, 'extrinsicProp');
                exprops = params(exprop_filter);
                if ~isempty(exprops)
                    if opts.FilterDirectDescendentExProps
                        parent_filter = arrayfun(@(p) ismember(p.Parent, C), exprops); % Ensure we only combine extrinsic props from direct descendents.  
                        sys_exprops = Combine(exprops(parent_filter)); % Combine extrinsic props from components into system level extrinsic props
                    else
                        sys_exprops = Combine(exprops);
                    end
                    for i = 1:numel(sys_exprops)
                        sys_exprops(i).Parent = obj;
                    end
                    params = [unique(params(~exprop_filter)); existing_params; unique(exprops); sys_exprops]; % Order is [component ~extrinsic props, system props, component extrinsic props, system extrinsic props]
                else
                    params = unique(params); % Select params with unique Sym_ID identifiers
                end
            end
            obj.Params = params;
        end
        
        function comp_arrays = Replicate(obj_array, N, opts)
            % Replicate is used to create N independent copies of all the components in obj_array
            % Returns comp_arrays cell array containing the duplicate components.  
            % - comp_arrays{i} is an 1xN Component array of copies corresponding to obj_array(i)
            arguments
                obj_array
                N {mustBeInteger}
                opts.ReplicateMethod ComponentReplicateMethods = "Copy"
                opts.Rename = true;
                opts.RedefineParams = false;
                opts.RedefineElement = true;
                opts.RedefineChildren = true;
            end
               
            if opts.ReplicateMethod == "Copy"
                comp_arrays = cell(size(obj_array));
                for i = 1:numel(obj_array)
                    comp = obj_array(i);
                    new_comps = Component.empty(N,0);
                    for j = 1:N
                        new_comp = copy(comp);
                        if opts.Rename
                            new_comp.Name = comp.Name+"_"+string(j);
                        end
                        if opts.RedefineParams
                            new_comp.DefineParams();
                            obj.init(); % Concrete method that can be overriden by subclasses
                        end
                        if opts.RedefineElement
                            new_comp.DefineElement();
                        end
                        if opts.RedefineChildren
                            new_comp.DefineChildren();
                        end
                        
                        new_comps(j,1) = new_comp;
                    end
                    comp_arrays{i} = new_comps;
                end
            elseif opts.ReplicateMethod == "Construct"
                comp_arrays = cell(size(obj_array));
                for i = 1:numel(obj_array)
                    comp = obj_array(i);
                    unique_props = setdiff(properties(comp), properties(class(comp))); % Get properties unique to the specific component, i.e. remove properties in the Component superclass
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

        function x = eq(obj1,obj2)
            % eq() method to compare a hetergeneous
            % array of SystemElements
            % GraphVertices are 'handle' objects, so they are considered
            % equivalent if the handles point to the same object.
            if numel(obj1)==numel(obj2)
                x = false(size(obj1));
                for i = 1:numel(obj1)
                    x(i) = eq@handle(obj1(i), obj2(i));
                end
            elseif numel(obj1)==1
                x = false(size(obj2));
                for i = 1:numel(obj2)
                    x(i) = eq@handle(obj1, obj2(i));
                end
            elseif numel(obj2)==1
                x = false(size(obj1));
                for i = 1:numel(obj1)
                    x(i) = eq@handle(obj2, obj1(i));
                end
            else
                error('array sizes are incompatible') 
            end
        end
        
        function x = ne(obj1, obj2)
            x = ~eq(obj1, obj2);
        end
    end
    
    %% Methods to be Overwritten by Subclasses
    methods
        function init_super(obj)
            DefineParams(obj);
            init(obj); % Optionally defined by subclasses
            DefineElement(obj);
            DefineChildren(obj);
        end
        
        function DefineParams(obj)
            cp = compParam.gatherObjectParams(obj);

            for i = 1:numel(cp)
                cp(i).Parent = obj;
            end
            
            obj.Params = cp;
        end
        
        function init(obj)
        end
        
        function DefineElement(obj)
        end
        
        function DefineChildren(obj)
        end
    end
end

