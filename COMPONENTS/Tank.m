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
        
        % Block Name
        Name char ='Tank'
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
            C(1) = Type_Capacitance("10"); 
            C(2) = Type_Capacitance("1"); 
             
            % Power Flow Types
            P(1) = Type_PowerFlow('Val_Char',"c*u*xt");
            P(2) = Type_PowerFlow('Val_Char',"c*(u1-u2)*xt");
            P(3) = Type_PowerFlow('Val_Char',"(u1-u2)");
            
            % define vertices
            Vertex(1) = GraphVertex_Internal('Description','Liquid Temp','Type',1,'Capacitance',C(1));
            Vertex(2) = GraphVertex_Internal('Description','Mass','Type',1,'Capacitance',C(2));
            Vertex(3) = GraphVertex_External('Description','Inlet');
            Vertex(4) = GraphVertex_External('Description','Outlet');
            Vertex(5) = GraphVertex_External('Description','Sink');
            
            % define edges
            Edge(1) = GraphEdge_Internal('PowerFlow',P(1),'Input',1,'Port',1);
            Edge(2) = GraphEdge_Internal('PowerFlow',P(1),'Input',2,'Port',2);
            Edge(3) = GraphEdge_Internal('PowerFlow',P(2),'Input',[1 2]);
            Edge(4) = GraphEdge_Internal('PowerFlow',P(3),'Input',[1 2]);
            Edge(5) = GraphEdge_External();
             
            g = Graph(E,Vertex,Edge);
            
        end
    end
end
