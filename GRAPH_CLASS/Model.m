classdef Model < matlab.mixin.Copyable
    %TYPE Summary of this class goes here
    %   Detailed explanation goes here
   
    
    properties 
        % System representation
        % x_dot = f(x,u,d)
        % y     = g(x,u,d)
        f_sym (:,1) sym
        g_sym (:,1) sym
        
        
        
    end
    
    properties (SetAccess = protected)
        
        LinearModel LinearizedModel = LinearizedModel.empty() % this could be another object

        
        CalcF (1,1) function_handle = @(x)0 % calculates x_dot
        CalcG (1,1) function_handle = @(x)0% calculates y
        
        Nx % number of states
        Nu % number of inputs
        Nd % number of disturbances
        Ny % number of outputs
    end
    
  
    
    methods
        function obj = Model(obj)
            a = 1;
        end

        function init(obj)
            a = 1;
            
        end
        
        function LinearizeModel(obj)
            
            x1       = sym('x%d'    ,[obj.Nx    1]); % dynamic states
            u1       = sym('u%d'    ,[obj.Nu    1]); % inputs
            d1       = sym('d%d'    ,[obj.Nd    1]);
            
            A = jacobian(obj.f_sym,x1);
            B = jacobian(obj.f_sym,u1);
            G = jacobian(obj.f_sym,d1);
            
            C = jacobian(obj.g_sym,x1);
            D = jacobian(obj.g_sym,u1);
            H = jacobian(obj.g_sym,d1);
            
            obj.LinearModel = LinearizedModel(A,B,G,C,D,H);
    
        end
        


        function set.f_sym(obj,val)
            obj.f_sym = val;
            updateCalcF(obj);
        end
        function updateCalcF(obj)          
            x1       = sym('x%d'    ,[obj.Nx    1]); % dynamic states
            u1       = sym('u%d'    ,[obj.Nu    1]); % inputs
            d1       = sym('d%d'    ,[obj.Nd    1]);
            
            obj.CalcF = matlabFunction(obj.f_sym,'Vars',[{[x1] [u1], [d1]}]);        
        end
        
        function set.g_sym(obj,val)
            obj.g_sym = val;
            updateCalcG(obj);
        end
        function updateCalcG(obj)          
            x1       = sym('x%d'    ,[obj.Nx    1]); % dynamic states
            u1       = sym('u%d'    ,[obj.Nu    1]); % inputs
            d1       = sym('d%d'    ,[obj.Nd    1]);
            
            obj.CalcG = matlabFunction(obj.g_sym,'Vars',[{[x1] [u1], [d1]}]);        
        end
        
      
        
    end
end

