classdef ComponentData
    %Value class encapsulating data corresponding to a physical component
    % It contains an array of compParamValues in its Data property that are
    % used to replace the values of system element parameters.  
    
    properties
        Component string % Corresponds to Component class name
        Make string
        Model string
        SerialNumber string
        Description string
        
        Data (:,1) compParamValue
    end
    
    methods
        function obj = ComponentData(comp, make, model, data)
            obj.Component = comp;
            obj.Make = make;
            obj.Model = model;
            obj.Data = data(:);
            
            for i = 1:numel(obj.Data)
                obj.Data(i).Component = obj.Component;
            end
        end
    end
end

