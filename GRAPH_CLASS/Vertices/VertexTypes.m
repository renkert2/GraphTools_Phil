classdef VertexTypes
    % VertexTypes is an enumeration class that combines Variable Type and 
    % Domain information in the Graph Modeling Toolbox. VertexTypes are 
    % used to check compatiabilty between component ports. 
    % Combines Domain and VariableType information into single enumeration
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Contributors: Christopher T. Aksland and Phil Renkert
    % Association: University of Illionis at Urbana-Champaign
    % Contact: aksland2@illinois.edu and renkert2@illinois.edu
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % potential improvements:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    enumeration
        % Domains derived from: https://www.mathworks.com/help/physmod/simscape/ug/basic-principles-of-modeling-physical-networks.html#bq89sba-3
        Abstract ('Abstract', 'Abstract')
        
        % Electrical Domain
        Voltage ('Electrical','Effort')
        Current ('Electrical','Flow')
        
        % Hydraulic Domain
        GaugePressure ('Hydraulic','Effort')
        VolumetricFlowRate ('Hydraulic','Flow')
        
        % Isothermal Liquid
        AbsolutePressure ('IsothermalLiquid','Effort')
        MassFlowRate ('IsothermalLiquid','Flow')
        
        % Magnetic Domain
        MagnetomotiveForce ('Magnetic','Effort')
        Flux ('Magnetic','Flow')
        
        % Mechanical Rotational Domain
        AngularVelocity ('MechanicalRotational','Effort')
        Torque ('MechanicalRotational','Flow')
        
        % Mechanical Translational Domain
        TranslationalVelocity ('MechanicalTranslational','Effort')
        Force ('MechanicalTranslational','Flow')
        
        % Thermal Domain
        Temperature ('Thermal','Effort')
        HeatFlow ('Thermal','Flow')
    end
    
    properties
        Domain Domains = 'Abstract'
        VariableType VariableTypes = 'Abstract'
    end
    
    methods
        function obj = VertexTypes(domain, var_type)
            obj.Domain = domain;
            obj.VariableType = var_type;
        end
        
        function x = isAbstract(obj)
            x = obj.VariableType == VariableTypes.Abstract;
        end
        
        function x = isEffort(obj)
            x = obj.VariableType == VariableTypes.Effort;
        end
        
        function x = isFlow(obj)
            x = obj.VariableType == VariableTypes.Flow;
        end
        
        function x = isConnectable(varargin) 
            % isConnectable = true -> Vertices can be connceted via an edge
            % Arguments are two VertexType arguments as comma separated list or array of two objects
            
            if nargin > 1
                obj_array = [varargin{:}];
            else
                obj_array = varargin{1};
            end
            
            assert(numel(obj_array) == 2, 'isConnectable requires two VertexType arguments as comma separated list or array of two objects');
            obj1 = obj_array(1);
            obj2 = obj_array(2);
            
            if (obj1.VariableType == VariableTypes.Abstract || obj2.VariableType == VariableTypes.Abstract)
                x = true;
            elseif (obj1.VariableType ~= obj2.VariableType)
                x = true;
            else
                x = false;
            end
        end
        
        function x = isCompatible(varargin) 
            % isCompatible = true -> Vertices can be combined in equivalence connection
             % Arguments are VertexTypes as comma separated list or array of VertexType objects
             
            if nargin > 1
                obj_array = [varargin{:}];
            else
                obj_array = varargin{1};
            end
                
            assert(numel(obj_array)>=2,'Array of two or more objects required');
            
            definedTypes = obj_array(VertexTypes.Abstract ~= obj_array); % Allow combination of abstract types
            if ~isempty(definedTypes)
                x = all(definedTypes(1) == definedTypes);
            else
                x = true;
            end
        end            
    end
        
end

