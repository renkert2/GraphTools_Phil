classdef Type_PowerFlow < Type
    %Type used to store power flow information
    
    methods
        function obj = Type_PowerFlow(type)
            [vars, params] = calcVars(type);
            obj = obj@Type(vars, params, type);
        end
        
        function updateInputs(obj, type)
            [vars, params] = calcVars(type);
            
            obj.vars = vars;
            obj.params = params;
        end
        
        function init(obj, type)
            obj.updateInputs(type)
            init@Type(obj, type);
        end
    end
            
end


function num_inputs = parseInputVars(type)
if isa(type, 'sym') || isa(type, 'char')
    type = string(type);
end

assert(isa(type, 'string'), "PowerFlow Type must be of type string, char, or sym");

patt = '(?<=u)\d+';

input_list = double(regexp(type, patt, 'match'));

num_inputs = max(input_list);
end

function [vars, params] = calcVars(type)
state_vars = [sym('xt'), sym('xh')];
num_inputs = parseInputVars(type);
if num_inputs > 0
    input_vars = sym('u', [1 num_inputs]);
    vars = [state_vars, input_vars];
    params = {state_vars(1), state_vars(2), input_vars};
else
    vars = [state_vars];
    params = {state_vars(1), state_vars(2)};
end

end

