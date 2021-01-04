% params
cp_f1 = 1;
cp_f2 = 1;
HTC = 1;


% edge matrix
E = [3 1; ...
    1 5; ...
    4 2; ...
    2 6; ...
    1 2];

% Capacitance Types
C(1) = Type_Capacitance("1");

% Power Flow Types
P(1) = Type_PowerFlow("u1*xt");
P(2) = Type_PowerFlow("xt");
P(3) = Type_PowerFlow("xh");

% Inputs
I(1) = GraphInput('Description', 'Input 1');
I(2) = GraphInput('Description', 'Input 2');

% Define Vertices
Vertex(1) = GraphVertex_Internal('Description','Temp S1','Type',1,'Capacitance',C(1),'Coefficient',10,'Initial',25);
Vertex(2) = GraphVertex_Internal('Description','Temp S2','Type',1,'Capacitance',C(1),'Coefficient',10,'Initial',25);
Vertex(3) = GraphVertex_External('Description','Inlet S1','Capacitance',C(1));
Vertex(4) = GraphVertex_External('Description','Inlet S2','Capacitance',C(1));
Vertex(5) = GraphVertex_External('Description','Outlet S1','Capacitance',C(1));
Vertex(6) = GraphVertex_External('Description','Outlet S2','Capacitance',C(1));

% Define Edges
Edge(1) = GraphEdge_Internal('PowerFlow',P(1),'Input',I(1),'Port',1,'Coefficient',cp_f1,'TailVertex',Vertex(E(1,1)),'HeadVertex',Vertex(E(1,2)));
Edge(2) = GraphEdge_Internal('PowerFlow',P(1),'Port',2,'Coefficient',cp_f1,'TailVertex',Vertex(E(2,1)),'HeadVertex',Vertex(E(2,2)));
Edge(3) = GraphEdge_Internal('PowerFlow',P(1),'Input',I(2),'Port',3,'Coefficient',cp_f2,'TailVertex',Vertex(E(3,1)),'HeadVertex',Vertex(E(3,2)));
Edge(4) = GraphEdge_Internal('PowerFlow',P(1),'Input',I(2),'Port',4,'Coefficient',cp_f2,'TailVertex',Vertex(E(4,1)),'HeadVertex',Vertex(E(4,2)));
Edge(5) = GraphEdge_Internal('PowerFlow',[P(2); P(3)], 'Input', [I(1); I(2)], 'Coefficient',[HTC -HTC],'TailVertex',Vertex(E(5,1)),'HeadVertex',Vertex(E(5,2)));

g = Graph(Vertex,Edge);

%% 

gModel = GraphModel(g)
gModel.init()