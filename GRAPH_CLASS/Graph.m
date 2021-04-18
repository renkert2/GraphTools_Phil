classdef Graph < matlab.mixin.Copyable
    % Graph is a class that stores all the information required to generate
    % a graph model. Graphs can be instatiated as an empty object or
    %
    % g = Graph(Vertex,Edge) where Vertex is an array of GraphVertex 
    % elements and Edge is array of GraphEdge elements.
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Contributors: Christopher T. Aksland and Phil Renkert
    % Association: University of Illionis at Urbana-Champaign
    % Contact: aksland2@illinois.edu and renkert2@illinois.edu
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % potential improvements:
    % - Consider added Matlab "digraph" or adjcancy matrix as a propertey.
    % - Consider improving how we find the head and tail vertex indices.
    % - At some point, graph can be defined with only edges
    % - Move the Internal/External etc get functions into the GraphVertex and
    %   GraphEdge classes. Since those functions only operate on vertices
    %   and edges, it seems they fit better there.
    % - The graph combination code should be moved to a separate file
    %   (files) since it's so long and complex.
    % - Phil: Figure out a more elegant solution for SymParams
    % - Use VPA on all our symbolic calculations to save a ton of time
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    properties
        Vertices (:,1) GraphVertex % set of vertices
        Edges (:,1) GraphEdge % set of edges
        
        Outputs (:,1) GraphOutput % Graph Outputs
    end
    
    properties (SetAccess = private)
        Inputs (:,1) GraphInput % graph Inputs
        
        Nv (1,1) double % Number of Internal Vertices
        Nx (1,1) double % Number of Dynamic States / Dynanic Vertices
        Ne (1,1) double % Number of Edges
        Nu (1,1) double % Number of Inputs 
        Nev (1,1) double % Number of External Vertices
        Nee (1,1) double % Number of External Edges
          
        E (:,2) double % Graph Edge Matrix
        M (:,:) double % Incidence Matrix
    end
    
    properties (Dependent)  
        v_tot (1,1) double % total number of vertices
  
        Tails (1,1) double 
        Heads (1,1) double 
    end
    
    properties (Dependent)
        InternalVertices  (:,1) GraphVertex
        DynamicVertices   (:,1) GraphVertex
        AlgebraicVertices (:,1) GraphVertex
        ExternalVertices  (:,1) GraphVertex
        InternalEdges     (:,1) GraphEdge
        ExternalEdges     (:,1) GraphEdge 
    end
    
    properties (SetAccess = ?SystemElement)
        Parent = Component.empty();
    end
    
    methods
        function obj = Graph(varargin) % constructor method
            % this function can be initialized with an edge matrix.
            % Immediately calculate the incidence matrix
           if nargin == 0
               % do nothing
           elseif nargin == 2
               obj.Vertices = varargin{1};
               obj.Edges = varargin{2};
               obj.init()
           else
              error('A Graph object must be initialized as an empty object or as Graph(Vertex_Set, Edge_Set ).') 
           end
                              
        end
                      
        function x = get.v_tot(obj) % total number of vertices
            x = obj.Nv + obj.Nev;
        end
        
        function x = get.Tails(obj) % tail incidence matrix
            x = double(obj.M'== 1);
        end
        
        function x = get.Heads(obj) % head incidence matrix
            x = double(obj.M'== -1);
        end
        
        function x = get.InternalVertices(obj) % get list of Internal Vertices
            x = getInternal(obj.Vertices);
        end
        
        function x = get.ExternalVertices(obj) % get list of external vertices
            x = getExternal(obj.Vertices);
        end 
        
        function x = get.DynamicVertices(obj) % get list of Dynamic Vertices
            x = getDynamic(obj.Vertices);
        end
        
        function x = get.AlgebraicVertices(obj) % get list of algebraic Vertices
            x = getAlgebraic(obj.Vertices);
        end
        
        function x = get.InternalEdges(obj) % get list of internal edges
            x = getInternal(obj.Edges);
        end
        function x = get.ExternalEdges(obj) % get list of external edges
            x = getExternal(obj.Edges);
        end
        
        function init(obj)
            obj.Vertices = [obj.DynamicVertices; obj.AlgebraicVertices; obj.ExternalVertices];
            obj.Edges    = [obj.InternalEdges; obj.ExternalEdges];
   
            % Calculate input array from Internal Edges
            obj.Inputs = unique(vertcat(obj.InternalEdges.Input), 'stable');

            % calculate graph size, order, etc
            obj.Nv  = length(obj.InternalVertices);
            obj.Nx  = length(obj.DynamicVertices);
            obj.Ne  = length(obj.InternalEdges);
            obj.Nu  = length(obj.Inputs);
            obj.Nev = length(obj.Vertices)-obj.Nv;
            obj.Nee = length(obj.Edges)-obj.Ne;
            
            % list of tail and head vertices
            vTail = vertcat(obj.InternalEdges.TailVertex);
            vHead = vertcat(obj.InternalEdges.HeadVertex);
            
            % indicies of tail and head vertices
            vT_idx =  arrayfun(@(x) find(arrayfun(@(y) eq(x,y), obj.Vertices)),vTail);
            vH_idx =  arrayfun(@(x) find(arrayfun(@(y) eq(x,y), obj.Vertices)),vHead);

            % calculate Edge and Incidence matrix
            obj.E = [vT_idx vH_idx];
            m = zeros(max(max(obj.E)),size(obj.E,1));
            for i = 1:size(obj.E,1)
                m(obj.E(i,1),i) =  1; % tails
                m(obj.E(i,2),i) = -1; % heads
            end
            obj.M = m;
            
        end
               
        function h = plot(obj,varargin)
            % basic digraph plotting.
            E_Augmented = obj.E; % edge matrix
            try % try to augment the edge matrix to include external edges
                E_idx = arrayfun(@(x) find(x==obj.InternalVertices),vertcat(obj.ExternalEdges.HeadVertex)); % index of vertices with external edges
                Eext = [[obj.v_tot+1:1:obj.v_tot+obj.Nee]' E_idx]; % edge matrix of external edges
                E_Augmented = [E_Augmented; Eext]; % augment E matrix with external edges
                skipPlotExt = 0;
            catch
                skipPlotExt = 1;
            end
            G = digraph(E_Augmented(:,1),E_Augmented(:,2)); % make digraph
            h = plot(G,varargin{:}); % plot digraph
            labeledge(h,E_Augmented(:,1)',E_Augmented(:,2)',[1:obj.Ne, 1:obj.Nee]); % edge lables
            highlight(h,[obj.Nv+1:1:obj.v_tot],'NodeColor','w') % remove external vertices
            xLoc = h.XData(obj.Nv+1:1:obj.v_tot); yLoc = h.YData(obj.Nv+1:1:obj.v_tot); % get external vertex location
            hold on; scatter(xLoc,yLoc,5*h.MarkerSize,'MarkerEdgeColor',h.NodeColor(1,:)); hold off; % add external vertices back
            if ~skipPlotExt
                highlight(h,reshape(Eext',1,[]),'LineStyle','--') % make external edges dashed
                highlight(h,[obj.v_tot+1:1:obj.v_tot+length(E_idx)],'NodeLabelColor','w','NodeColor','w') % remove augmented external edge tail vertices and labels
                xLoc = [ h.XData(Eext(:,1))]; yLoc = [ h.YData(Eext(:,1))]; % get external edge tail vertex location
            end
            hold on; scatter(xLoc,yLoc,5*h.MarkerSize,'MarkerEdgeColor',h.EdgeColor(1,:)); hold off; % plot external edge tail vertices
        end
        
        function [obj, ConnectE, ConnectV] = Combine(G, ConnectE, ConnectV, opts) % Input is vector of GraphClass elements
            % G - GraphClass Class Array
            % ConnectX - Cell vector with first list containing component
            % indices; second list containing equivalent properties (edges
            % or vertices) corresponding to the components.
            
            arguments
                G (:,1) Graph
                ConnectE cell = {}
                ConnectV cell = {}
                opts.CopyEdges logical = true;
            end
            
            % Format ConnectV and ConnectE as {:,1} cell array containing list of
            % equivalent vertices and edges, respectively
            
            if ~isempty(ConnectE)
                ConnectE = formatConnectX(ConnectE, 'GraphEdge', 'Edges');
                Nce = length(ConnectE); % Number of edge connections
                Neconn_delta = cellfun(@length, ConnectE);% Number of edges involved in each connection
                Neconn = sum(Neconn_delta); % Total edges involved in connections
                
                ConnectV_E = cell(2*Nce, 1);
                edge_conn_map(Nce,2) = GraphEdge();
            else
                Nce = 0;
                Neconn_delta = 0;
                Neconn = 0;
                
                ConnectV_E = {};
                edge_conn_map = GraphEdge.empty();
            end
            
            if nargin == 3 % If ConnectV is specified
                ConnectV = formatConnectX(ConnectV, 'GraphVertex', 'Vertices');
            end
            
            input_conn_map = GraphInput.empty();
            
            %% Parse ConnectE for necessary Vertex Connections, create edge_conn_map, and create input_conn_map:
            % - Appends ConnectV with necessary vertex connections resulting from edge connections
            % edge_conn_map:
            % - First column contains edges replaced in connection that will be modified in the connection,
            % - Second column contains primary edges resulting from connection
            % - s.t. edge_conn_map(:,1) becomes edge_conn_map(:,2)
            % input_conn_map:
            % - First column contains inputs replaced in connection that will be modified in the connection,
            % - Second column contains primary inputs resulting from connection
            % - s.t. input_conn_map(:,1) becomes input_conn_map(:,2)
           
            for ce = 1:Nce
                edges = fliplr(ConnectE{ce});
                edge_conn_map(ce, :) = edges;
                                   
                equiv_heads = [edges.HeadVertex];
                assert(isCompatible([equiv_heads.VertexType]), 'Incompatible head vertices in edge connection %d', ce) % Check for compatible head vertices with isCompatible method of VertexTypes
                
                equiv_tails = [edges.TailVertex];
                assert(isCompatible([equiv_tails.VertexType]), 'Incompatible tail vertices in edge connection %d', ce) % Check for compatible tail vertices with isCompatible method of VertexTypes
                
                ConnectV_E{2*ce - 1} = equiv_heads;
                ConnectV_E{2*ce} = equiv_tails;
                
                old_inputs = [edges(1).Input];
                new_inputs = [edges(2).Input];
                
                if (~isempty(old_inputs)) && (~isempty(new_inputs))
                    if numel(old_inputs) == numel(new_inputs)
                        for i = 1:numel(old_inputs)
                            input_conn_map(end+1,:) = [old_inputs(i) new_inputs(i)]; % Input from first edge in ConnectE will replace input from second edge in ConnectE
                        end
                    else
                        error('Incompatible inputs in edge connection %d', ce);
                    end
                end
            end
            
            if isempty(ConnectV)
                ConnectV = ConnectV_E;
            else
                ConnectV = [ConnectV; ConnectV_E];
            end
            
            %% Create vertex_conn_map:
            % - First column contains all vertices that will be modified in the connection,
            % - Second column contains primary vertices that verts in the first column will map to
            % - s.t. Vertex_Conn_Map(:,1) becomes Vertex_Conn_Map(:,2)
            
            % Determine the Primary Vertex for each connection
            Ncv = length(ConnectV); % Number of vertex connections
            Nvconn_delta = cellfun(@length, ConnectV); % Number of vertices involved in each connection
            Nvconn = sum(Nvconn_delta); % Total vertices involved in connections
            
            primary_vertices(Ncv,1) = GraphVertex();
            vertex_conn_map(Nvconn-Ncv, 2) = GraphVertex();
            vconnmap_counter = cumsum([0; Nvconn_delta-1]);
            for cv = 1:length(ConnectV)
                verts = ConnectV{cv};
                
                % Check for compatible Vertex Types with isCompatible method of VertexTypes
                types = [verts.VertexType];
                assert(isCompatible(types), 'Incompatible VertexTypes in Vertex Connection %d of ConnectV', cv);
                
                % Check Number of Internal Vertices and Assign Primary Vertex
                int_vert_flags = arrayfun(@(x) isa(x,'GraphVertex_Internal'), verts);
                n_int_verts = sum(int_vert_flags);
                if n_int_verts == 0
                    primary_vertex = verts(1);
                elseif n_int_verts == 1
                    primary_vertex = verts(int_vert_flags);
                else
                    error('Error in Vertex Connection %d.  Only one Internal Vertex can be involved in a connection.', cv)
                end
                primary_vertices(cv,1) = primary_vertex;
                
                conn_verts = verts(primary_vertex ~= verts);
                
                r = vconnmap_counter(cv)+(1:numel(conn_verts)); % Range of elements in vertex_conn_map to be assigned
                vertex_conn_map(r, 1) = conn_verts;
                vertex_conn_map(r,2) = primary_vertex;
            end
            
            %% Process vertex_conn_map
            vertex_conn_map = unique(vertex_conn_map, 'rows','stable'); % remove duplicate mappings in vertex_conn_map
            if length(vertex_conn_map(:,1)) ~= length(unique(vertex_conn_map(:,1)))
                error('Vertices assigned to multiple primary vertices.');
            end
            
            % Reformat serial connections in vertex_conn_map
            vertex_conn_map = ReformatSerialConnections(vertex_conn_map);
            
            %% Process input_conn_map
            if ~isempty(input_conn_map)
                input_conn_map = unique(input_conn_map, 'rows','stable'); % remove duplicate mappings in input_conn_map
                if length(input_conn_map(:,1)) ~= length(unique(input_conn_map(:,1)))
                    error('Inputs assigned to different primary inputs.');
                end
                
                % Reformat serial connections in input_conn_map
                input_conn_map = ReformatSerialConnections(input_conn_map);
            end
            %% Construct Vertex and Edge Vectors
            all_verts = vertcat(G.Vertices);
            all_edges = vertcat(G.Edges);
            all_outputs = vertcat(G.Outputs);
           
            sys_verts = all_verts(~ismember(all_verts, vertex_conn_map(:,1)));
            
            if ~isempty(edge_conn_map)
                sys_edges = all_edges(~ismember(all_edges, edge_conn_map(:,1)));
            else
                sys_edges = all_edges;
            end
            
            if opts.CopyEdges
                sys_edges = copy(sys_edges);
            end
            
            for i = 1:numel(sys_edges)
                % Replace head vertex in system edges
                head_v = sys_edges(i).HeadVertex;
                
                log_index = head_v == vertex_conn_map(:,1);
                if any(log_index)
                    primary_vertex = vertex_conn_map(log_index, 2);
                    sys_edges(i).HeadVertex = primary_vertex;
                end
                

                if isa(sys_edges(i), 'GraphEdge_Internal')
                    % Replace tail vertex in system edges
                    tail_v = sys_edges(i).TailVertex;
                    
                    log_index = tail_v == vertex_conn_map(:,1);
                    if any(log_index)
                        primary_vertex = vertex_conn_map(log_index, 2);
                        sys_edges(i).TailVertex = primary_vertex;
                    end
                    
                    % Replace inputs in system edges
                    if ~isempty(input_conn_map)
                        inputs = [sys_edges(i).Input];
                        if ~isempty(inputs)
                            for j = 1:numel(inputs)
                                log_index = inputs(j) == input_conn_map(:,1);
                                if any(log_index)
                                    sys_edges(i).Input(j) = input_conn_map(log_index,2);
                                end
                            end
                        end
                    end
                end
            end
            
            for i = 1:numel(all_outputs)                    
                    % Replace inputs in system outputs
                    if ~isempty(input_conn_map)
                        inputs = [];
                        idx = [];
                        for j = 1:length(all_outputs(i).Breakpoints)
                            if isa(all_outputs(i).Breakpoints{j},'GraphInput')
                                inputs = [inputs all_outputs(i).Breakpoints{j}];
                                idx = [ idx j];
                            end
                        end
                        if ~isempty(inputs)
                            for j = 1:numel(inputs)
                                log_index = inputs(j) == input_conn_map(:,1);
                                if any(log_index)
                                    all_outputs(i).Breakpoints{idx(j)} = input_conn_map(log_index,2);
                                end
                            end
                        end
                    end
            end
            
            
            % Instantiate and Return the Graph Object
            obj = Graph(sys_verts, sys_edges);
            obj.Outputs = all_outputs;
         
            %% Helper Functions
            function ConnectX = formatConnectX(ConnectX, class, prop)
                if all(cellfun(@(x) isa(x, class), ConnectX),'all') % Connect X Already given as lists of equivalent GraphVertices or GraphEdges with dominant element in front of list
                    return
                elseif all(cellfun(@(x) isa(x, 'double'), ConnectX),'all') % ConnectX specified with component indices in 1st row, Property indices in second row
                    num_connections = size(ConnectX,2);
                    if size(ConnectX,1) == 3 % If third row specifying dominant component is given
                        ConnectX_temp = cell(2, num_connections);
                        for c = 1:num_connections
                            comps = ConnectX{1,c};
                            props = ConnectX{2,c};
                            dom_comp = ConnectX{3,c};
                            
                            i = dom_comp == comps ; % Identify dominant component in first row
                            
                            new_comps = [comps(i) comps(~i)];
                            new_props = [props(i) props(~i)];
                            
                            ConnectX_temp{1,c} = new_comps;
                            ConnectX_temp{2,c} = new_props;
                        end
                        ConnectX = ConnectX_temp;
                    end
                    
                    vec_lengths = cellfun(@numel, ConnectX);
                    assert(all(vec_lengths(1,:)==vec_lengths(2,:)), 'Component and Property Vectors must be the same length for each connection');
                    
                    ConnectX_temp = cell(num_connections, 1);
                    for c = 1:num_connections
                        graphs_i = ConnectX{1,c};
                        elems_i = ConnectX{2,c};
                        elems = arrayfun(@(c,e) G(c).(prop)(e), graphs_i, elems_i, 'UniformOutput', false); % Arrayfun can't handle heterogenous object arrays so nonuniform output required
                        elems = [elems{:}]; % Convert cell array to heterogenous object array
                        ConnectX_temp{c} = elems;
                    end
                    ConnectX = ConnectX_temp;
                else
                    error('Invalid ConnectX cell array');
                end
            end
            
            function conn_map = ReformatSerialConnections(conn_map)
                inters = intersect(conn_map(:,1),conn_map(:,2)); % Find common elements in first and second columns of conn_map
                if ~isempty(inters) % If intersections exist
                    for i = 1:numel(inters)
                        target_val = conn_map(inters(i) == conn_map(:,1),2); 
                        conn_map(inters(i) == conn_map(:,2),2) = target_val;
                    end
                end
            end
        end
    end
end

