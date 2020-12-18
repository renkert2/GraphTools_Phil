classdef Type_PowerFlow < Type
    %Type used to store power flow information
    
    methods
        function obj = Type_PowerFlow(varargin)
            obj = obj@Type(varargin{:});
            obj.init()
        end
        
        function init(obj)
            
            % list of all symbolic variables
            T_var_all = symvar(obj.Val_Sym).';
            T_var_red = T_var_all;
            
            % get list of number of feasible inputs (Options) and defined
            % inputs (T_var_red). inputs are special here because an 
            % arbitrary number of inputs can be used
            T_var_red(ismember(T_var_red,[sym('xh');sym('xt')])) = []; % remove head and tail states from the list of options
            Options = sym('u',[length(T_var_red) 1]); % total number of allowable inputs
            
            % check to see if T_var is a subset of variable options
            if sum(ismember(T_var_red,Options)) == length(T_var_red)
                obj.vars = T_var_all;
                params = [{sym('xt') sym('xh') sym('u',[1 max(length(T_var_red),2)])}]; %[xt xh u] % max function necessary to make sure the the input parameter can be variable size
                obj.Val_Func = matlabFunction(obj.Val_Sym,'Vars',params);
                obj.Jac_Sym  = jacobian(obj.Val_Sym,[sym('xt') sym('xh') sym('u',[1 length(T_var_red)])]);
                obj.Jac_Func = matlabFunction(obj.Jac_Sym,'Vars',params);
            else
                error('Invalid powerflow definition. Define power flows in term of tail state xt, head state xh, and inputs u1, u2, ... uN')
            end
            
%             init@Type(obj);
        end
    end
end
    

