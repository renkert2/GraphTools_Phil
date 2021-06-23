function y = Model_SimulinkInterpretedFunction(u, sys_name, obj_name, func)
%MODEL_SIMULINKINTERPRETEDFUNCTION Summary of this function goes here
%   Detailed explanation goes here
    mdlWks = get_param(sys_name, 'ModelWorkspace');
    m = evalin(mdlWks,obj_name); % Import Model object from ModelWorkspace
    y = func(m, u);
end