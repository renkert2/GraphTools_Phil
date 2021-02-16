classdef LinearModel < matlab.mixin.Copyable
    %LINEARMODEL Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        % System Description
        % x_dot = A*x + B*u + G*d + f0
        % y     = C*x + D*u + H*d + g0
        
        A_sym  (:,:) sym
        B_sym  (:,:) sym
        E_sym  (:,:) sym
%         f0_sym (:,:) double % these come from the Model class for now...
        
        C_sym  (:,:) sym
        D_sym  (:,:) sym
        H_sym  (:,:) sym
%         g0_sym (:,:) double % these come from the Model class for now...
    end
    
    properties
        % System Description
        % x_dot = A*x + B*u + G*d + f0
        % y     = C*x + D*u + H*d + g0
        
        CalcState  (1,1) function_handle = @(x)0
        CalcOutput  (1,1) function_handle = @(x)0
%         CalcE  (1,1) function_handle = @(x)0
%         CalcF0 (:,:) double % these come from the Model class for now...
        
%         CalcC  (1,1) function_handle = @(x)0
%         CalcD  (1,1) function_handle = @(x)0
%         CalcH  (1,1) function_handle = @(x)0
%         CalcG0 (:,:) double % these come from the Model class for now...
    end
    
    methods
    
        function obj = LinearModel(varargin)
            if nargin == 6
%                 x1       = sym('x%d'    ,[obj.Nx    1]); % dynamic states
%                 u1       = sym('u%d'    ,[obj.Nu    1]); % inputs
%                 d1       = sym('d%d'    ,[obj.Nd    1]);
                obj.A_sym = varargin{1};
                obj.B_sym = varargin{2};
                obj.E_sym = varargin{3};
                obj.C_sym = varargin{4};
                obj.D_sym = varargin{5};
                obj.H_sym = varargin{6};
                
                
            else
                
            end
        end
        
        
    end
end

