classdef Model < matlab.mixin.Copyable
    %TYPE Summary of this class goes here
    %   Detailed explanation goes here

    properties (SetAccess = protected) 
        Nx % number of states
        Nu % number of inputs
        Nd % number of disturbances      
        Ny % number of outputs

    end
    
    properties 
        % System representation
        % x_dot = f(x,u,d)
        % y     = g(x,u,d)
        f_sym (:,1) sym
        g_sym (:,1) sym
          
    end
    
    properties (Dependent)
        StateNames (:,1) string
        InputNames (:,1) string
        DisturbanceNames (:,1) string
        OutputNames (:,1) string    
    end
    
    properties (SetAccess = protected)
        LinModel LinearModel = LinearModel.empty() % this could be another object
        
        CalcF_func (1,1) function_handle = @(x)0 % calculates x_dot
        CalcG_func (1,1) function_handle = @(x)0% calculates y  
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
            %x1       = sym('x%d'    ,[obj.Nx    1]); % dynamic states
            x1 = genSymVars('x%d', obj.Nx);
            %u1       = sym('u%d'    ,[obj.Nu    1]); % inputs
            u1 = genSymVars('u%d', obj.Nu);
            %d1       = sym('d%d'    ,[obj.Nd    1]); % disturbances
            d1 = genSymVars('d%d', obj.Nd);
            
            A = jacobian(obj.f_sym,x1);
            B = jacobian(obj.f_sym,u1);
            E = jacobian(obj.f_sym,d1);
            
            C = jacobian(obj.g_sym,x1);
            D = jacobian(obj.g_sym,u1);
            H = jacobian(obj.g_sym,d1);
            
            obj.LinModel = LinearModel(A,B,E,C,D,H);
            
            obj.CalcF_func = matlabFunction(obj.f_sym,'Vars',[{[x1], [u1], [d1]}]);
            obj.CalcG_func = matlabFunction(obj.g_sym,'Vars',[{[x1], [u1], [d1]}]);

            obj.LinModel.CalcState  = matlabFunction(obj.LinModel.A_sym,obj.LinModel.B_sym,obj.LinModel.E_sym,'Vars',[{[x1], [u1], [d1]}]);
            obj.LinModel.CalcOutput = matlabFunction(obj.LinModel.C_sym,obj.LinModel.D_sym,obj.LinModel.H_sym,'Vars',[{[x1], [u1], [d1]}]);
%             obj.LinearModel.CalcOutput = matlabFunction(obj.LinearModel.B_sym,'Vars',[{[x1] [u1], [d1]}]);
%             obj.LinearModel.CalcE = matlabFunction(obj.LinearModel.E_sym,'Vars',[{[x1] [u1], [d1]}]);
%             obj.LinearModel.CalcC = matlabFunction(obj.LinearModel.C_sym,'Vars',[{[x1] [u1], [d1]}]);
%             obj.LinearModel.CalcD = matlabFunction(obj.LinearModel.D_sym,'Vars',[{[x1] [u1], [d1]}]);
%             obj.LinearModel.CalcH = matlabFunction(obj.LinearModel.H_sym,'Vars',[{[x1] [u1], [d1]}]);
                     
        end
        
        function F = CalcF(obj,x,u,d)
            assert((size(x,2) == size(u,2)) && (size(x,2) == size(d,2)),"x, u, and d require equivalent number of columns")
            if size(x,2) == 1
                F = obj.CalcF_func(x,u,d);
            else
                F = splitapply(obj.CalcF_func,x,u,d,1:size(x,2));
            end
        end
         
        function G = CalcG(obj,x,u,d)
            assert((size(x,2) == size(u,2)) && (size(x,2) == size(d,2)),"x, u, and d require equivalent number of columns")
            if size(x,2) == 1
                G = obj.CalcG_func(x,u,d);
            else
                G = splitapply(obj.CalcG_func,x,u,d,1:size(x,2));
            end
        end
        
        function [A,B,E,F0,C,D,H,G0] = Linearize(obj,x0,u0,d0)
            
            % 
            [A,B,E] = obj.LinModel.CalcState(x0,u0,d0);
%             B  = obj.LinearModel.CalcB(x0,u0,d0);
%             E  = obj.LinearModel.CalcE(x0,u0,d0);
            F0 = obj.CalcF_func(x0,u0,d0) - A*x0 - B*u0 - E*d0;
            
            [C,D,E]  = obj.LinModel.CalcOutput(x0,u0,d0);
%             D  = obj.LinearModel.CalcD(x0,u0,d0);
%             G  = obj.LinearModel.CalcG_func(x0,u0,d0);
            G0 = obj.CalcG_func(x0,u0,d0) - C*x0 - D*u0 - G*d0;
               
        end
        
        function x = get.StateNames(obj)
            x = defineStateNames(obj);
        end
        
        function x = get.InputNames(obj)
            x = defineInputNames(obj);
        end
        
        function x = get.DisturbanceNames(obj)
            x = defineDisturbanceNames(obj);                       
        end
        
        function x = get.OutputNames(obj)
            x = defineOutputNames(obj);                       
        end   
    end
    
    methods (Abstract)
        defineStateNames
        defineInputNames
        defineDisturbanceNames
        defineOutputNames
    end
end

