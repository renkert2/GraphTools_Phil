classdef Tank < Component
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
        function obj = Tank(varargin)          
            obj@Component(varargin{:}); % calls the superclass constructor                    
        end
    end
    
    methods (Access = protected)
        function DefineComponent(obj)
            
            % edge matrix
            E = [3 1; ...
                 1 4; ...
                 1 5; ...
                 5 2];
             
            % Capacitance Types
            C(1) = Type_Capacitance("1");
             
            % Power Flow Types
            P(1) = Type_PowerFlow("u1*xt"); 
            P(2) = Type_PowerFlow("(u1-u2)*xt");
            P(3) = Type_PowerFlow("(u1-u2)");
          
            % define vertices
            Vertex(1) = GraphVertex_Internal('Description','Liquid Temp','Capacitance',C(1),'Coefficient',obj.cp_f,'Initial',25, 'VertexType', 'Temperature');
            Vertex(2) = GraphVertex_Internal('Description','Tank Mass','Capacitance',C(1),'Coefficient',1,'Initial',100);
            Vertex(3) = GraphVertex_External('Description','Inlet');
            Vertex(4) = GraphVertex_External('Description','Outlet');
            Vertex(5) = GraphVertex_External('Description','Sink');
%             Vertex(3) = GraphVertex_External('Description','Inlet','Capacitance',C(1));
%             Vertex(4) = GraphVertex_External('Description','Outlet','Capacitance',C(1));
%             Vertex(5) = GraphVertex_External('Description','Sink','Capacitance',C(1));            
            % define Inputs
            I(1) = GraphInput('Description',"Inflow",'Bounds',Limits(0,1));
            I(2) = GraphInput('Description',"Outflow",'Bounds',Limits(0,1));
            
            % define edges
            Edge(1) = GraphEdge_Internal('PowerFlow',P(1),'Input',I(1),'Coefficient',obj.cp_f,'TailVertex',Vertex(E(1,1)),'HeadVertex',Vertex(E(1,2)));
            Edge(2) = GraphEdge_Internal('PowerFlow',P(1),'Input',I(2),'Coefficient',obj.cp_f,'TailVertex',Vertex(E(2,1)),'HeadVertex',Vertex(E(2,2)));
            Edge(3) = GraphEdge_Internal('PowerFlow',P(2),'Input',[I(1) I(2)],'Coefficient',obj.cp_f,'TailVertex',Vertex(E(3,1)),'HeadVertex',Vertex(E(3,2)));
            Edge(4) = GraphEdge_Internal('PowerFlow',P(3),'Input',[I(1) I(2)],'Coefficient',1,'TailVertex',Vertex(E(4,1)),'HeadVertex',Vertex(E(4,2)));
            Edge(5) = GraphEdge_External('HeadVertex',Vertex(1),'Description','Heat Load');
            
            % special lookup functions
            T = Type('x1'); 
            Vertex(1).CapFunction = LookupFunction('Function',T,'Breakpoints',[Vertex(2)]);
            
            % Build Graph
            g = Graph(Vertex,Edge);
            obj.Graph = g;
            
            % Define Ports
            p(1) = ComponentPort('Description','Inflow','Element',Edge(1));
            % skip port 2 since it's a source in the simulink block
            p(3) = ComponentPort('Description','Outflow','Element',Edge(2));
            obj.Ports = p;
            
            % Define Additional Outputs
            O = Type('b1*b2+b3', {sym('b1') sym('b2') sym('b3')});
            o = GraphOutput('Description','Test','Function',O,'Breakpoints',{Vertex(1) I(2) Vertex(2)});
            obj.Graph.Outputs = o;
            
        end
    end
end
