classdef GraphVertex_Internal < GraphVertex
    % GraphVertex_Internal defines the properties of internal vertices in a
    % graph object in the Graph Modeling Toolbox
    % Instatiate an empty object or use an input parser
    % Definable properties (with class) include:
    % - Description (string)
    % - VertexType (VertexTypes)
    % - DynamicType (DynamicTypes)
    % - CapFunction (LookupFunction)
    % - Capacitance (Type_Capacitance) 
    % - Coefficient ({mustBeNonnegativeOrSym})
    % - Initial (double)
    % - Bounds (Limits) 
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Contributors: Christopher T. Aksland and Phil Renkert
    % Association: University of Illionis at Urbana-Champaign
    % Contact: aksland2@illinois.edu and renkert2@illinois.edu
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % potential improvements:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    properties
        CapFunction (:,1) LookupFunction
        Capacitance (:,1) Type_Capacitance
        Coefficient (:,1) {mustBeNonnegativeOrSym} = 0
        Initial (1,1) double = 0
        Bounds Limits {mustBeScalarOrEmpty}
    end
end

