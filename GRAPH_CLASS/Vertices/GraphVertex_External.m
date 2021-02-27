classdef GraphVertex_External < GraphVertex
    % GraphVertex_External defines the properties of external vertices in a
    % graph object in the Graph Modeling Toolbox
    % Instatiate an empty object or use an input parser
    % Definable properties (with class) include:
    % - Description (string)
    % - VertexType (VertexTypes)
    % - DynamicType (DynamicTypes)
    % - (no additional)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Contributors: Christopher T. Aksland and Phil Renkert
    % Association: University of Illionis at Urbana-Champaign
    % Contact: aksland2@illinois.edu and renkert2@illinois.edu
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % potential improvements:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    methods
        function obj = GraphVertex_External(varargin) % constructor method
            obj@GraphVertex(varargin{:}); % calls the superclass constructor
        end
         
    end

end

