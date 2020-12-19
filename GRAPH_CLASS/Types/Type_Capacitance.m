classdef Type_Capacitance < Type
    %Type used to store vertex capacitance information
    methods
        function obj = Type_Capacitance(type)
            vars = [sym('x')];
            params = {vars}; 
            obj = obj@Type(vars, params, type);
        end
    end
end

