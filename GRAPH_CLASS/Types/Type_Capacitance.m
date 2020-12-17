classdef Type_Capacitance < Type
    %Type used to store vertex capacitance information
    
    methods
        function obj = Type_Capacitance(varargin)
            obj = obj@Type(varargin{:});
            obj.init()
        end
        
        function init(obj)
            % list of all symbolic variables
            T_var_all = symvar(obj.Val_Sym).';
%             T_var_red = T_var_all;
            
            % get list of number of feasible states (Options) and defined
            % inputs (T_var_red).
%             T_var_red(ismember(T_var_red,[sym('xh');sym('xt')])) = []; % remove head and tail states from the list of options
            Options = sym('x'); % total number of allowable inputs
            
            % check to see if T_var is a subset of variable options
            if sum(ismember(T_var_all,Options)) == length(T_var_all)
                obj.vars = T_var_all;
                params = [{sym('x')}]; %[x]
                obj.Val_Func = matlabFunction(obj.Val_Sym,'Vars',params);
                obj.Jac_Sym  = jacobian(obj.Val_Sym,[sym('x')]);
                obj.Jac_Func = matlabFunction(obj.Jac_Sym,'Vars',params);
            else
                error('Invalid capacitance definition. Define capacitance in term of state x.')
            end
            
        end
    end
end
    
