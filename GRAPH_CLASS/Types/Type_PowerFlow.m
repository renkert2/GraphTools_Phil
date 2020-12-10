classdef Type_PowerFlow < Type
    %Type used to store power flow information
    
    methods
        function obj = Type_PowerFlow(varargin)
            obj = obj@Type(varargin{:});
            obj.init()
        end
        
        function init(obj)
            % use symvar(obj.Val_Sym) to get most of these
            % at this point, maybe have error checking that power flows
            % were defined in terms of valid variables xt, xh, u1, ..., uN
            syms x_t x_h u
            syms c g 
            syms dg_dxt dg_dxh dg_du
            obj.vars = [x_t x_h u c g dg_dxt dg_dxh];
%             init@Type(obj);
        end
    end
end
    

