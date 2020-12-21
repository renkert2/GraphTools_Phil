classdef LinearizedModel < Model %matlab.mixin.Copyable
    %LINEARMODEL Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        % System Description
        % x_dot = A*x + B*u + G*d + f0
        % y     = C*x + D*u + H*d + g0
        
        A_sym  (:,:) double
        B_sym  (:,:) double
        G_sym  (:,:) double
%         f0_sym (:,:) double % these come from the Model class for now...
        
        C_sym  (:,:) double
        D_sym  (:,:) double
        H_sym  (:,:) double
%         g0_sym (:,:) double % these come from the Model class for now...
    end
    
    properties
        % System Description
        % x_dot = A*x + B*u + G*d + f0
        % y     = C*x + D*u + H*d + g0
        
        CalcA  (:,:) double
        CalcB  (:,:) double
        CalcG1  (:,:) double
%         CalcF0 (:,:) double % these come from the Model class for now...
        
        CalcC  (:,:) double
        CalcD  (:,:) double
        CalcH  (:,:) double
%         CalcG0 (:,:) double % these come from the Model class for now...
    end
    
    methods
    
        function obj = LinearizedModel(varargin)
            if nargin == 6
                x1       = sym('x%d'    ,[obj.Nx    1]); % dynamic states
                u1       = sym('u%d'    ,[obj.Nu    1]); % inputs
                d1       = sym('d%d'    ,[obj.Nd    1]);
                
                
                
            else
                
            end
        end
        
        
    end
end

