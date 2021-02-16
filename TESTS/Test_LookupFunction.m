clear all
close all

C(1) = Type_Capacitance("1");
Vertex(1) = GraphVertex_Internal('Description','Liquid Temp','Capacitance',C(1),'Coefficient',[1000 .5],'Initial',25, 'VertexType', 'Temperature');
Vertex(2) = GraphVertex_Internal('Description','Tank Mass','Capacitance',C(1),'Coefficient',1,'Initial',100);
I(1) = GraphInput("Tank Input 1");



T = Type([sym('p1') sym('p2')],{sym('p1') sym('p2')},'p1*p2');
LF = LookupFunction('Function',T,'Breakpoints',{Vertex(2) I(1)});
Vertex(1).CapFunction = LF;


