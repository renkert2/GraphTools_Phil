function y = Model_SimulinkInterpretedFunction(u,model_name,func)
%MODEL_SIMULINKINTERPRETEDFUNCTION Summary of this function goes here
%   Detailed explanation goes here
    mdlWks = get_param(model_name, 'ModelWorkspace');
    bm = evalin(mdlWks,'ModelObject');
    y = func(bm, u);
end