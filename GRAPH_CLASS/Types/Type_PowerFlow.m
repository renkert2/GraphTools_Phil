classdef Type_PowerFlow < Type
    %Type used to store power flow information
    
    methods
        function obj = Type_PowerFlow(varargin)
            obj = obj@Type(varargin{:});
            obj.init()
        end
        
        function init(obj)
            syms x_t x_h u
            syms c g 
            syms dg_dxt dg_dxh dg_du
            obj.vars = [x_t x_h u c g dg_dxt dg_];
            init@Type(obj);
        end
    end
end
    

