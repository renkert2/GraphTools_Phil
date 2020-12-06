classdef HeatExchanger < Component_Super
    % HeatExchanger is a class the defines a heat exchanger model
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Author: Christopher T. Aksland
    % Association: University of Illionis at Urbana-Champaign
    % Contact: aksland2@illinois.edu
    % Revision History:
    % 7/6/2020 - Class creation
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    properties
        
        % Block Name
        Name char ='Heat Exchanger'
        % Side 1 working Fluid
        fluid1 char = 'JP8'
        % Side 2 working Fluid
        fluid2 char = 'water'
        % Initial Side 1 Fluid temperature [C]
        T1_init(1,1) double {mustBeNumeric} = 25;
        % Initial Side 2 Fluid temperature [C]
        T2_init(1,1) double {mustBeNumeric} = 25;
        % Side 1 fluid Specific Heat [J/kg]
        cp_f1 (1,1) double {mustBeNumeric} = 2000;
        % Side 2 fluid Specific Heat [J/kg]
        cp_f2 (1,1) double {mustBeNumeric} = 2000;
        % Heat Transfer Coefficient [W/K]
        HTC (1,1) double {mustBeNumeric} = 10;
        
    end
    
    methods
        function obj = HeatLoad(varargin)
            
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

