classdef GraphOutput < handle
    % GraphOutput is a class that describes additional outputs for a graph
    % in the Graph Modeling Toolbox. GraphOutputs can be instatiated as an 
    % empty object or use the input parser.
    % Definable properties (with class) include:
    % - Description (string)
    % - Function (Type or symfun)
    % - Breakpoints (cell array of GraphVertices or GraphInputs)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Contributors: Christopher T. Aksland and Phil Renkert
    % Association: University of Illionis at Urbana-Champaign
    % Contact: aksland2@illinois.edu and renkert2@illinois.edu
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % potential improvements:
    % Creat a breakpoints class that superclasses GraphVertices and Graph
    %   Inputs
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    properties
        Description (1,1) string = "Default"
        Function % Must be symfun or Type, asserted in set method
        Breakpoints cell % Must be cell array of GraphVertices or GraphInputs, asserted in set method
    end
    
    properties (SetAccess = ?Component)
        Parent Component = Component.empty()
    end
    
    methods
        function obj = GraphOutput(varargin) % constructor method
            if nargin > 1
                obj = my_inputparser(obj,varargin{:});
            end
        end
        
        function set.Function(obj, val) % error checking on function property
            assert(isa(val, 'Type') || isa(val, 'symfun'), 'Function must be a Type object or symfun');
            obj.Function = val;
        end
        
        function set.Breakpoints(obj,val) % error checking on Breakpoints property
            for i = 1:length(val)
                if ~isa(val{i},'GraphVertex') && ~isa(val{i},'GraphInput')
                    error('Breakpoints must be of type GraphVertex or GraphInput')
                end
            end
            obj.Breakpoints = val;
        end
        
    end
end

