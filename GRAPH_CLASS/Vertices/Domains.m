classdef Domains
    enumeration
        Abstract
        Electrical
        Hydraulic
        IsothermalLiquid
        Magnetic
        MechanicalRotational
        MechanicalTranslational
        Thermal
    end
    
    methods
        function x = isCompatible(varargin)
            if nargin > 1
                obj_array = [varargin{:}];
            else
                obj_array = varargin{1};
            end
            
            assert(numel(obj_array)>=2,'Array of two or more Domain objects required');
            
            definedDomains = obj_array(Domains.Abstract ~= obj_array); % Allow combination of abstract domains
            if ~isempty(definedDomains)
                x = all(definedDomains(1) == definedDomains);
            else
                x = true;
            end
        end
    end
end