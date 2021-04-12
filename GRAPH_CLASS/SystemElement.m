classdef SystemElement < matlab.mixin.Heterogeneous & handle
    %SYSTEMELEMENT Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Name (1,1) string = "Component" % Block Name
        Graph (1,1) Graph
        Ports (:,1) ComponentPort = ComponentPort.empty()
        
        extrinsicProps (:,1) extrinsicProp
        SymParams SymParams {mustBeScalarOrEmpty}
    end
    
    methods
        function obj = SystemElement(inputArg1,inputArg2)
            %SYSTEMELEMENT Construct an instance of this class
            %   Detailed explanation goes here
            obj.Property1 = inputArg1 + inputArg2;
        end
        
        function outputArg = method1(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
    end
end

