function [Sys, PropMaps] = GenSysGraph_Phil(G, ConnectV, ConnectE) % Input is vector of GraphClass elements
% G - GraphClass Class Array
% ConnectX - Cell vector with first list containing component
% indices; second list containing equivalent properties (edges
% or vertices) corresponding to the components.

assert(all(class(G)=='Graph'), 'Input G must be GraphModel array');
assert(iscolumn(G), 'Input G must be column array');

% Vertex Property Map
[V_map, Nv, Nt] = VertexPropMap();
Nx = Nv-Nt;

% Edge Property Map
[E_map,Ne] = EdgePropMap();

Vertices(Nv) = GraphVertex;
Edges(Ne) = GraphEdge;

for g_i = length(G):-1:1 % Loop backwards so earlier components overide later ones
    Vertices(V_map{g_i}) = G(g_i).Vertices;
    
    Edges_temp = G(g_i).Edges;
    Edges_temp.HeadVertex = 
    Edges_temp.TailVertex = 
    Edges(E_map{g_i}) = G(g_i).Edges;
end

PropMaps = struct();
PropMaps.V_map = V_map;
PropMaps.E_map = E_map;
PropMaps.U_map = U_map;

% Set System Graph Properties
Sys = GraphModel(Ne,Nv,Nt,Nu);
fn = fieldnames(temp_g);
for i = 1:length(fn)
    Sys.(fn{i}) = temp_g.(fn{i});
end

% Call the GraphModel init function to calculate F_inner,etc
Sys.init();

    function [V_map, Nv, Nt] = VertexPropMap()
        sink_flags = cellfun(@(v_conn) all(v_conn{2}>arrayfun(@(i) G(i).Nx, v_conn{1})), ConnectV); % Connections resulting in sink vertices
        
        [chi_Conn, chi_Conn_sink] = ConnectedProps(ConnectV, sink_flags);
        assert(~CountRepeats(chi_Conn), 'ConnectV cannot have subvertices used in multiple vertex equivalencies.  Please write as one connection');
        
        % Vertex Property Filters: filters to tell us which partition of system vertex vector
        % chi_s a vertex belongs to
        V_filters = cell(size(G)); % Tells us which partition each vertex belongs to for each component
        N_partitions = 4; % Number of partitions of resulting system vertex vector chi_s
        for g_i = 1:length(G)
            v = 1:G(g_i).Nv;
            
            % Find vertices that are being connected
            connection_filter = ismember(v,chi_Conn{g_i});
            
            % State filter: find vertices that are state vertices
            state_filter = v<=G(g_i).Nx;
            
            %res sink filter - find what vertices map to sink state vertices
            res_sink_filter = ismember(v,chi_Conn_sink{g_i});
            
            filter = zeros([N_partitions G(g_i).Nv]);
            %chi_s_ubar_filter - State vertices not involved in connections
            filter(1,:) = state_filter & ~connection_filter;
            %chi_s_hat_ubar_filter - Vertices connected to create state vertices
            filter(2,:) = ~res_sink_filter & connection_filter;
            %chi_s_hat_lbar_filter - Vertices connected to create sink vertices
            filter(3,:) = res_sink_filter & connection_filter;
            %chi_s_lbar_filter = Sink vertices not involved in connections
            filter(4,:) = ~state_filter & ~connection_filter;
            
            V_filters{g_i} = filter;
        end
        
        N_chi_s_ubar = sum(cellfun(@(x) sum(x(1,:)), V_filters));%sum([G.Nx]'-chi_Conn_ubar); % Calculates the total number of state vertices not from connections
        N_chi_s_lbar = sum(cellfun(@(x) sum(x(4,:)), V_filters)); % Calculates the total number of sink vertices not from connections
        
        N_chi_s_hat = length(ConnectV);% - CountRepeats(chi_Conn);
        N_chi_s_hat_lbar = sum(sink_flags);% - CountRepeats(chi_Conn_sink); % Calculates total number of sink vertices from connetions
        N_chi_s_hat_ubar = N_chi_s_hat-N_chi_s_hat_lbar; % Calculates the total number of state vertices from connections
        
        N_chi_s = [N_chi_s_ubar, N_chi_s_hat_ubar, N_chi_s_hat_lbar, N_chi_s_lbar];
        
        Nv = sum(N_chi_s);
        Nt = N_chi_s_hat_lbar + N_chi_s_lbar;
        
        % Vertex Property Map: Rows = Graphs, Columns = vertex index that corresponding
        % subsystem vertex maps to vertex in the resulting system
        partition_counter =  zeros(1,4);% Initialize counter to keep track of the number of elements in each partition
        part_index = cumsum([0 N_chi_s(1:end-1)]);
        V_map = cell(size(G));
        for g_i = 1:length(G)
            V_map{g_i} = zeros(1,G(g_i).Nv);
        end
        for g_i = 1:length(G)
            for p_i = [1,4] % Loop over chi_s partitions 1 and 4 not involving connections
                filter = V_filters{g_i}(p_i,:);
                indices = filter.*cumsum(filter);
                V_map{g_i} = V_map{g_i}+(part_index(p_i)+partition_counter(p_i)).*filter+indices;
                partition_counter(p_i) = partition_counter(p_i)+sum(filter);
            end
        end
        for c_i = 1:length(ConnectV) % Loop over connections (chi_s partitions 2 and 3)
            if sink_flags(c_i)
                p_i = 3;
            else
                p_i=2;
            end
            v_i = part_index(p_i)+partition_counter(p_i)+1;
            comps = ConnectV{c_i}{1};
            verts = ConnectV{c_i}{2};
            for i = 1:length(comps)
                V_map{comps(i)}(verts(i)) = v_i;
            end
            partition_counter(p_i) = partition_counter(p_i)+1;
        end
    end
    function [E_map, Ne] = EdgePropMap()
        N_e= [G.Ne]';
        Ne = sum(N_e)-length(ConnectE); % Number of resulting system edges
        N_xi_s_hat = length(ConnectE); % Each edge connection results in a single edge
        N_xi_s = Ne - N_xi_s_hat; % Number of edges not resulting from connection
        
        [chi_Conn, ~] = ConnectedProps(ConnectE, zeros(1,length(ConnectE)));
        assert(~CountRepeats(chi_Conn), 'The same edge cannot be used in multiple edge equivalencies');
        
        E_filters = cell(size(G)); % Tells us which partition each vertex belongs to for each component
        N_partitions = 2; % Number of partitions of resulting system vertex vector chi_s
        for g_i = 1:length(G)
            e = 1:N_e(g_i);
            
            % Find edges that are being connected
            connection_filter = ismember(e,chi_Conn{g_i});
            
            filter = zeros([N_partitions N_e(g_i)]);
            %chi_s_ubar_filter - edges not involved in connections
            filter(1,:) = ~connection_filter;
            %chi_s_hat_ubar_filter - edges involved in connections
            filter(2,:) = connection_filter;
            
            E_filters{g_i} = filter;
        end
        
        partition_counter = zeros(1,N_partitions);
        part_index = [0 N_xi_s];
        E_map = cell(size(G));
        for g_i = 1:length(G)
            E_map{g_i} = zeros(1,N_e(g_i));
        end
        p_i = 1; % xi_s partition 1, not involving connections
        for g_i = 1:length(G)
            filter = E_filters{g_i}(p_i,:); % Indicates which component edges are not involved in connections
            indices = filter.*cumsum(filter);
            E_map{g_i} = E_map{g_i}+(part_index(p_i)+partition_counter(p_i)).*filter+indices;
            partition_counter(p_i) = partition_counter(p_i)+sum(filter);
        end
        p_i = 2;  % xi_s partition 2, involving connections
        for c_i = 1:length(ConnectE)
            e_i = part_index(p_i)+partition_counter(p_i)+1;
            comps = ConnectE{c_i}{1};
            edges = ConnectE{c_i}{2};
            for i = 1:length(comps)
                E_map{comps(i)}(edges(i)) = e_i;
            end
            partition_counter(p_i) = partition_counter(p_i)+1;
        end
    end
    function [U_map, Nu] = InputPropMap() % Consider equivalent inputs - use edge map to see if they're equivalent.  Important for fluid flow advection
        U_map = cell(size(G));
        u_cnt = 0;
        for g_i=1:length(G)
            U_map{g_i} = (1:G(g_i).Nu)+u_cnt;
            u_cnt = u_cnt+G(g_i).Nu;
        end
        Nu = u_cnt;
    end
    function [ConnectedProps, ConnectedPropsFiltered] = ConnectedProps(ConnectX, filter)
        % Transforms ConnectX input into a cell vector where each
        % element is a vector of properties (vertices or edges)
        % involved in a ConnectX connection.  Optional filter is a
        % list of flags corresponding to ConnectX.
        if nargin==1
            filter = ones(size(G));
        end
        
        ConnectedProps = cell(size(G));
        ConnectedPropsFiltered = cell(size(G));
        
        for g_i = 1:length(G)
            x = [];
            x_filtered = [];
            for i = 1:length(ConnectX)
                comps = ConnectX{i}{1};
                Xs = ConnectX{i}{2};
                if ismember(g_i, comps)
                    x = [x Xs(comps==g_i)];
                    if filter(i)
                        x_filtered = [x_filtered Xs(comps==g_i)];
                    end
                end
            end
            ConnectedProps{g_i}=x;
            ConnectedPropsFiltered{g_i}=x_filtered;
        end
    end
    function n_repeat = CountRepeats(chiConn)
        accum = arrayfun(@(x) accumarray(x,1), chiConn, 'UniformOutput',false);
        accum = vertcat(accum{:});
        n_repeat = sum(accum(accum>1)-1);
    end
end