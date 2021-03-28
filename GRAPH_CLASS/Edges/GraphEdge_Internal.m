classdef GraphEdge_Internal < GraphEdge
    % GraphEdge_Internal defines the properties of internal edges in a
    % graph object in the Graph Modeling Toolbox
    % Instatiate an empty object or use an input parser
    % Definable properties (with class) include:
    % - Description (string)
    % - PowerFlow (Type_PowerFlow) 
    % - Input (GraphInput)
    % - Coefficient (double)
    % - TailVertex (GraphVertex) 
    % - HeadVertex (GraphVertex) 
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Contributors: Christopher T. Aksland and Phil Renkert
    % Association: University of Illionis at Urbana-Champaign
    % Contact: aksland2@illinois.edu and renkert2@illinois.edu
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % potential improvements:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    properties
        PowerFlow (:,1) Type_PowerFlow % PowerFlow representation for an edge
        Input (:,1) GraphInput % inputs incident on the edge 
        Coefficient (:,1) {mustBeNumericOrSym} = 0 % the constant edge coefficient.  Must be a positive double() or a sym() or symParam()
        TailVertex (:,1) GraphVertex {mustBeScalarOrEmpty} % the edge tail vertex
        HeadVertex (:,1) GraphVertex {mustBeScalarOrEmpty} % the edge head vertex
    end
end
