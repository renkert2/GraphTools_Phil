classdef Type_Capacitance < Type
    % Type_Capacitance (subclass of see Type) is used to define the 
    % expression to calculate constant coefficient capacitances in the 
    % Graph Modeling Toolbox. The capacitance must be defined in terms of 
    % the vertex state 'x'.
    % Instatiate an object in the following form:
    % 
    % T = Type_Capacitance(type) where type defines the expression used to
    % calulacte the capacitance of a vertex in a Graph object. type must
    % only be a function of the variable x (as a string or symbolic
    % variable).
    % Ex: Suppose the capacitance of a vertex is calculated as c1*(1+x+x^2)
    % where c1 is a constant coefficient and x is the state of the vertex.
    % String definition: 
    % T = Type_Capacitance('1+x+x^2');
    % Symbolic definition:
    % var = sym('x')
    % T = Type_Capacitance(1+var+var^2)
      
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Contributors: Christopher T. Aksland and Phil Renkert
    % Association: University of Illionis at Urbana-Champaign
    % Contact: aksland2@illinois.edu and renkert2@illinois.edu
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % potential improvements:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    methods
        function obj = Type_Capacitance(type) % object constructor
            params = {sym('x')}; 
            obj = obj@Type(type, params); % call parent constructor
        end
    end
end

