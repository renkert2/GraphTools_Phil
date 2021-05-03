classdef compParamValue
    % Value class used by compParam to load parameter values representing
    % physical components.  
    
    properties
        Sym string
        Value double
        Unit 
        Description string
        
        Component string % Name of component object corresponding to this parameter
        % May want to add ComponentData parent here later
    end
    
    methods
        function obj = compParamValue(sym, val, opts)
            arguments
                sym 
                val
                opts.Unit = ""
                opts.Description = ""
                opts.Component = ""
            end
            
            obj.Sym = sym;
            obj.Value = val;
            obj.Unit = opts.Unit;
            obj.Description = opts.Description;
            obj.Component = opts.Component;
        end
    end
end

