classdef HeatLoad < Component
    % HeatLoad is a class the defines a heat load model.
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Author: Christopher T. Aksland
    % Association: University of Illionis at Urbana-Champaign
    % Contact: aksland2@illinois.edu
    % Revision History:
    % 7/6/2020 - Class creation
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    properties

        % Working Fluid
        fluid char = 'JP8'
        % Initial Fluid temperature [C]
        T_init(1,1) double {mustBeNumeric} = 25;
        % Fluid Specific Heat [J/kg]
        cp_f (1,1) double {mustBeNumeric} = 2000;  
    end
    
    methods
        function obj = HeatLoad(varargin)          
            obj@Component(varargin{:}); % calls the superclass constructor                     
        end
    end
    
    methods (Access = protected)
        function DefineComponent(obj)
            % Edge Matrix
            E = [2 1; ...
                 1 3];
            
            % Capacitance Types
            C(1) = Type_Capacitance("1");
            
            % Power Flow Types
            P(1) = Type_PowerFlow("u1*xt");
            
            % Define Vertices
            Vertex(1) = GraphVertex_Internal('Description','Load Temp','Capacitance',C(1), 'VertexType', 'Temperature');
            Vertex(2) = GraphVertex_External('Description','Inlet');
            Vertex(3) = GraphVertex_External('Description','Outlet');
%                         Vertex(2) = GraphVertex_External('Description','Inlet','Capacitance',C(1));
%             Vertex(3) = GraphVertex_External('Description','Outlet','Capacitance',C(1));
            
            % Define Inputs
            I(1) = GraphInput('Description',"HeatLoad Input 1",'Bounds',Limits(0,1));
             
            % Define Edges
            Edge(1) = GraphEdge_Internal('PowerFlow',P(1),'Input',I(1),'Coefficient',obj.cp_f,'TailVertex',Vertex(E(1,1)),'HeadVertex',Vertex(E(1,2)));
            Edge(2) = GraphEdge_Internal('PowerFlow',P(1),'Input',I(1),'Coefficient',obj.cp_f,'TailVertex',Vertex(E(2,1)),'HeadVertex',Vertex(E(2,2)));
            Edge(3) = GraphEdge_External('HeadVertex',Vertex(1),'Description','Heat Load');

            % Build Graph
            g = Graph(Vertex,Edge);
            obj.Graph = g;
            
            % Define Ports
            p(1) = ComponentPort('Description','Inflow','Element',Edge(1));
            % skip port 2 since it's a source in the simulink block
            p(3) = ComponentPort('Description','Outflow','Element',Edge(2));
            obj.Ports = p;
            
            % Define Additional Outputs
            O = Type('b1*b2');
            o = GraphOutput('Description','Test','Function',O,'Breakpoints',{Vertex(1) I(1)});
            obj.Graph.Outputs = o;
            
        end
    end
end
