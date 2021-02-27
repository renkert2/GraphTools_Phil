classdef Domains
    % Domains is an enumeration class that defines different energy 
    % domains within the Graph Modeling Toolbox. Domains are used to check
    % compatiabilty between component ports. Supported domains are
    % - Abtract (general/default value)
    % - Electrical
    % - Hydraulic
    % - IsothermalLiquid
    % - Magnetic
    % - MechanicalRotational
    % - MechanicalTranslational
    % - Thermal
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Contributors: Christopher T. Aksland and Phil Renkert
    % Association: University of Illionis at Urbana-Champaign
    % Contact: aksland2@illinois.edu and renkert2@illinois.edu
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % potential improvements:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % @Phil can you add additional comments this since you're more familiar.
    
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
        function x = isCompatible(varargin) % checks compatiblity between two domains for interconnections
            if nargin > 1
                obj_array = [varargin{:}];
            else
                obj_array = varargin{1};
            end
            
            % throw a warning message incase only one domain is specified
            assert(numel(obj_array)>=2,'Array of two or more Domain objects required.');
            
            definedDomains = obj_array(Domains.Abstract ~= obj_array); % Allow combination of abstract domains
            if ~isempty(definedDomains)
                x = all(definedDomains(1) == definedDomains);
            else
                x = true;
            end
        end
    end
end