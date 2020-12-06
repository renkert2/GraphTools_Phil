classdef SplitJunction < Component_Super
    % SplitJunction is a class the defines a fluid split or junction model
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Author: Christopher T. Aksland
    % Association: University of Illionis at Urbana-Champaign
    % Contact: aksland2@illinois.edu
    % Revision History:
    % 7/6/2020 - Class creation
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    properties
        
        % Block Name
        Name char ='Split Junction'
        % Working Fluid
        fluid char = 'JP8'
        % Initial Fluid temperature [C]
        T_init(1,1) double {mustBeNumeric} = 25;
        % Fluid Specific Heat [J/kg]
        cp_f (1,1) double {mustBeNumeric} = 2000;
        % Number of inflows 
        n_in(1,1) double {mustBeInterger} = 1;
        % Number of outflows 
        n_out(1,1) double {mustBeInterger} = 1;
    end
    
    methods
        function obj = SplitJunction(varargin)
            
            % populate properties
            obj = my_inputparser(obj,varargin{:});
            
        end
        
        function obj = GraphModel(obj)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
    end
end

