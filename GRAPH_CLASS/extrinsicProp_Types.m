classdef extrinsicProp_Types
        
    enumeration
        Abstract (@sum, 1)
        Mass (@sum, 2)
        Cost (@sum, 3)
    end
    
    properties
        combFunc function_handle % function used to map enumeration values from Component -> System
        id uint8 % Used for functions like 'Sort'
    end
    
    methods
        function obj = extrinsicProp_Types(fun, id)
            obj.combFunc = fun;
            obj.id = id;
        end
        
        function sorted = sort(obj_array)
            [~,i] = sort([obj_array.id]);
            sorted = obj_array(i);
        end
    end

end
