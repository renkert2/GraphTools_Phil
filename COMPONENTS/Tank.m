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
        function g = DefineGraph(obj)
            
            % edge matrix
            E = [3 1; ...
                 1 4; ...
                 1 5; ...
                 5 2];
            % Capacitance Types
            C(1) = Type_Capacitance("1"); 
            C(2) = Type_Capacitance("x"); 
             
            % Power Flow Types
            P(1) = Type_PowerFlow("u1*xt");
            P(2) = Type_PowerFlow("(u1-u2)*xt");
            P(3) = Type_PowerFlow("(u1-u2)");
          
            % define vertices
            Vertex(1) = GraphVertex_Internal('Description','Liquid Temp','Type',1,'Capacitance',[C(2) C(1)],'Coefficient',[1000 .5],'Initial',25);
            Vertex(2) = GraphVertex_Internal('Description','Tank Mass','Type',1,'Capacitance',C(1),'Coefficient',1,'Initial',100);
            Vertex(3) = GraphVertex_External('Description','Inlet','Capacitance',C(1));
            Vertex(4) = GraphVertex_External('Description','Outlet','Capacitance',C(1));
            Vertex(5) = GraphVertex_External('Description','Sink','Capacitance',C(1));
            
            % define Inputs
            I(1) = GraphInput("Tank Input 1");
            I(2) = GraphInput("Tank Input 2");
            
            % define edges
            Edge(1) = GraphEdge_Internal('PowerFlow',P(1),'Input',I(1),'Port',1,'Coefficient',obj.cp_f,'TailVertex',Vertex(E(1,1)),'HeadVertex',Vertex(E(1,2)));
            Edge(2) = GraphEdge_Internal('PowerFlow',P(1),'Input',I(2),'Port',3,'Coefficient',obj.cp_f,'TailVertex',Vertex(E(2,1)),'HeadVertex',Vertex(E(2,2)));
            Edge(3) = GraphEdge_Internal('PowerFlow',P(2),'Input',[I(1) I(2)],'Coefficient',obj.cp_f,'TailVertex',Vertex(E(3,1)),'HeadVertex',Vertex(E(3,2)));
            Edge(4) = GraphEdge_Internal('PowerFlow',P(3),'Input',[I(1) I(2)],'Coefficient',1,'TailVertex',Vertex(E(4,1)),'HeadVertex',Vertex(E(4,2)));
            Edge(5) = GraphEdge_External('HeadVertex',Vertex(1),'Description','Heat Load');
             
            g = Graph(Vertex,Edge);
            
        end
    end
end
