classdef Type_Capacitance < Type
    %Type used to store vertex capacitance information
    
    methods
        function obj = Type_Capacitance(varargin)
            obj = obj@Type(varargin{:});
            obj.init()
        end
        
        function init(obj)
            syms x c g dg_dx
            obj.vars = [x c g dg_dx];
%             init@Type(obj);
        end
    end
end
    
