classdef LinearModel < matlab.mixin.Copyable
    % LinearModel is a propety of the model class and has linear state
    % space model matrices. The model is of the form:
    %
    % System Description
    % x_dot = A*x + B*u + G*d + f0
    % y     = C*x + D*u + H*d + g0
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Contributors: Christopher T. Aksland and Phil Renkert
    % Association: University of Illionis at Urbana-Champaign
    % Contact: aksland2@illinois.edu and renkert2@illinois.edu
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % potential improvements:
    % f0 and g0 should be moved into the linear model class...
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    properties      
        A_sym  (:,:) sym
        B_sym  (:,:) sym
        E_sym  (:,:) sym
        
        C_sym  (:,:) sym
        D_sym  (:,:) sym
        H_sym  (:,:) sym
    end
    
    properties
        % System Description
        % x_dot = A*x + B*u + G*d + f0
        % y     = C*x + D*u + H*d + g0
        
        CalcState function_handle {mustBeScalarOrEmpty} % outputs A, B, and G matrices
        CalcOutput function_handle {mustBeScalarOrEmpty} % outputs C, D, and H matrices
    end
    
    methods
        function obj = LinearModel(varargin) % constructor stores the linear matrix information
            if nargin == 6
                obj.A_sym = varargin{1};
                obj.B_sym = varargin{2};
                obj.E_sym = varargin{3};
                obj.C_sym = varargin{4};
                obj.D_sym = varargin{5};
                obj.H_sym = varargin{6};             
            end
        end  
    end
end

