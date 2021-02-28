classdef GraphInput < handle
    % GraphInput is a class that describes an input in the Graph Model
    % Toolbox. GraphInputs can be instatiated as an empty object, input
    % parsing, or as:
    % 
    % u = GraphInput(desc) where desc is the input description
    %
    % Definable properties (with class) include:
    % - Description (string)
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
        Description (1,1) string = "Default"
        Bounds Limits {mustBeScalarOrEmpty} 
    end
    
    properties (SetAccess = ?Component)
        Parent Component = Component.empty()
    end
    
    methods
        function obj = GraphInput(varargin)
            if nargin > 0
                if nargin == 1 && (isstring(varargin{1}) || ischar(varargin{1}))
                    desc = string(varargin{1});
                    obj.Description = desc;
                elseif nargin > 1
                    obj = my_inputparser(obj,varargin{:});
                else
                    error('Invalid argument to GraphInput constructor');
                end
            end
        end
    end
end

