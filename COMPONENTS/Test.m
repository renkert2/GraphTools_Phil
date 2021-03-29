classdef Test < Component
    % Tank is a class the defines a tank model
    
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
        % Initial tank mass [kg]
        m_init (1,1) double {mustBeNumeric} = 10;
        % Initial Tank temperature [C]
        T_init(1,1) double {mustBeNumeric} = 25;
        % Fluid Specific Heat [J/kg]
        cp_f (1,1) double {mustBeNumeric} = 2000;

    end
    
    methods
        function obj = Test(varargin)          
            obj@Component(varargin{:}); % calls the superclass constructor                    
        end
    end
    
    methods (Access = protected)
        function DefineComponent(obj)
            
            % edge matrix
            E = [5 1; 1 4; 6 2; 2 7; 1 3; 3 2; 4 8];
             
            % Capacitance Types
            C(1) = Type_Capacitance("1");
             
            % Power Flow Types
            P(1) = Type_PowerFlow("u1"); 
            P(2) = Type_PowerFlow("xt");
            P(3) = Type_PowerFlow("xh");
          
            % define vertices
            Vertex(1) = GraphVertex_Internal('Description','T1','Capacitance',C(1),'Coefficient',1,'Initial',100);
            Vertex(3) = GraphVertex_Internal('Description','W','Capacitance',C(1),'Coefficient',1,'Initial',100);
            Vertex(2) = GraphVertex_Internal('Description','T2','Capacitance',C(1),'Coefficient',1,'Initial',100);
            Vertex(4) = GraphVertex_Internal('Description','T12','Capacitance',C(1),'Coefficient',1,'Initial',100);
            Vertex(5) = GraphVertex_External('Description','Inlet 1');
            Vertex(8) = GraphVertex_External('Description','Outlet 1');
            Vertex(6) = GraphVertex_External('Description','Inlet 2');
            Vertex(7) = GraphVertex_External('Description','Outlet 2');
%             Vertex(3) = GraphVertex_External('Description','Inlet','Capacitance',C(1));
%             Vertex(4) = GraphVertex_External('Description','Outlet','Capacitance',C(1));
%             Vertex(5) = GraphVertex_External('Description','Sink','Capacitance',C(1));            
            % define Inputs
            I(1) = GraphInput('Description',"m1",'Bounds',Limits(0,1));
            I(2) = GraphInput('Description',"m2",'Bounds',Limits(0,1));
            
            % define edges
            Edge(1) = GraphEdge_Internal('PowerFlow',[P(1) P(2)],'Input',I(1),'Coefficient',[1 1],'TailVertex',Vertex(E(1,1)),'HeadVertex',Vertex(E(1,2)));
            Edge(2) = GraphEdge_Internal('PowerFlow',[P(1) P(2)],'Input',I(1),'Coefficient',[2 2],'TailVertex',Vertex(E(2,1)),'HeadVertex',Vertex(E(2,2)));
            Edge(3) = GraphEdge_Internal('PowerFlow',[P(1) P(2)],'Input',I(2),'Coefficient',[3 3],'TailVertex',Vertex(E(3,1)),'HeadVertex',Vertex(E(3,2)));
            Edge(4) = GraphEdge_Internal('PowerFlow',[P(1) P(2)],'Input',I(2),'Coefficient',[4 4],'TailVertex',Vertex(E(4,1)),'HeadVertex',Vertex(E(4,2)));
            Edge(5) = GraphEdge_Internal('PowerFlow',[P(2) P(3)],'Coefficient',[5 5],'TailVertex',Vertex(E(5,1)),'HeadVertex',Vertex(E(5,2)));
            Edge(6) = GraphEdge_Internal('PowerFlow',[P(2) P(3)],'Coefficient',[6 6],'TailVertex',Vertex(E(6,1)),'HeadVertex',Vertex(E(6,2)));
            Edge(7) = GraphEdge_Internal('PowerFlow',[P(1) P(2)],'Input',I(1),'Coefficient',[7 -7],'TailVertex',Vertex(E(7,1)),'HeadVertex',Vertex(E(7,2)));
                     
            % Build Graph
            g = Graph(Vertex,Edge);
            obj.graph = g;
            
%             % Define Ports
%             p(1) = ComponentPort('Description','Inflow','Element',Edge(1));
%             % skip port 2 since it's a source in the simulink block
%             p(3) = ComponentPort('Description','Outflow','Element',Edge(2));
%             obj.Ports = p;
            
        end
    end
end
