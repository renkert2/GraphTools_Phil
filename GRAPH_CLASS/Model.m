classdef Model < matlab.mixin.Copyable
    %TYPE Summary of this class goes here
    %   Detailed explanation goes here
   
    
        
    properties (SetAccess = protected)
        
        Nx % number of states
        Nu % number of inputs
        Nd % number of disturbances
%         Ny % number of outputs
    end
    
    properties 
        % System representation
        % x_dot = f(x,u,d)
        % y     = g(x,u,d)
        f_sym (:,1) sym
        g_sym (:,1) sym
          
    end
    
    properties (SetAccess = protected)
        
        LinModel LinearModel = LinearModel.empty() % this could be another object
        
        CalcF (1,1) function_handle = @(x)0 % calculates x_dot
        CalcG (1,1) function_handle = @(x)0% calculates y
        
    end
    
  
    
    methods
        function obj = Model(varargin)
            a = 1; 
            % add functionality here at some point
            
%             if nargin == 0
%                 % do nothing
%             elseif nargin == 1
%                 disp(sprintf('Model of class %s created.\n',class(varargin{1}))) % update this to list valid arguments
% 
%             elseif nargin == 2
%                 obj.f_sym = varargin{1};
%                 obj.g_sym = varargin{2};
%                 init(obj);
%             end
                
        end

        function init(obj)
            x1       = sym('x%d'    ,[obj.Nx    1]); % dynamic states
            u1       = sym('u%d'    ,[obj.Nu    1]); % inputs
            d1       = sym('d%d'    ,[obj.Nd    1]); % disturbances
            
            A = jacobian(obj.f_sym,x1);
            B = jacobian(obj.f_sym,u1);
            E = jacobian(obj.f_sym,d1);
            
            C = jacobian(obj.g_sym,x1);
            D = jacobian(obj.g_sym,u1);
            H = jacobian(obj.g_sym,d1);
            
            obj.LinModel = LinearModel(A,B,E,C,D,H);
            
            obj.CalcF = matlabFunction(obj.f_sym,'Vars',[{[x1], [u1], [d1]}]);
            obj.CalcG = matlabFunction(obj.g_sym,'Vars',[{[x1], [u1], [d1]}]);

            obj.LinModel.CalcState  = matlabFunction(obj.LinModel.A_sym,obj.LinModel.B_sym,obj.LinModel.E_sym,'Vars',[{[x1], [u1], [d1]}]);
            obj.LinModel.CalcOutput = matlabFunction(obj.LinModel.C_sym,obj.LinModel.D_sym,obj.LinModel.H_sym,'Vars',[{[x1], [u1], [d1]}]);
%             obj.LinearModel.CalcOutput = matlabFunction(obj.LinearModel.B_sym,'Vars',[{[x1] [u1], [d1]}]);
%             obj.LinearModel.CalcE = matlabFunction(obj.LinearModel.E_sym,'Vars',[{[x1] [u1], [d1]}]);
%             obj.LinearModel.CalcC = matlabFunction(obj.LinearModel.C_sym,'Vars',[{[x1] [u1], [d1]}]);
%             obj.LinearModel.CalcD = matlabFunction(obj.LinearModel.D_sym,'Vars',[{[x1] [u1], [d1]}]);
%             obj.LinearModel.CalcH = matlabFunction(obj.LinearModel.H_sym,'Vars',[{[x1] [u1], [d1]}]);
                     
        end
        
        function [A,B,E,F0,C,D,H,G0] = Linearize(obj,x0,u0,d0)
            
            % 
            [A,B,E] = obj.LinModel.CalcState(x0,u0,d0);
%             B  = obj.LinearModel.CalcB(x0,u0,d0);
%             E  = obj.LinearModel.CalcE(x0,u0,d0);
            F0 = obj.CalcF(x0,u0,d0) - A*x0 - B*u0 - E*d0;
            
            [C,D,E]  = obj.LinModel.CalcOutput(x0,u0,d0);
%             D  = obj.LinearModel.CalcD(x0,u0,d0);
%             G  = obj.LinearModel.CalcG(x0,u0,d0);
            G0 = obj.CalcG(x0,u0,d0) - C*x0 - D*u0 - G*d0;
               
        end
           
        
    end
end

