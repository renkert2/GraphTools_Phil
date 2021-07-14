classdef extrinsicProp < compParam
    % extrinsicProp is a class that represents design-dependent 
    % properties like mass and price.  When Components are 
    % combined into a system, extrinsicProps are combined into the
    % total system extrinsicProp.  E.X. if component one had a Mass extrinsicProp 
    % with value m1 and component 2 had a Mass extrinsicProp with value m2, 
    % the combined system would have a resulting Mass extrinsicProp with value
    % m1+m2
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Contributors: Christopher T. Aksland and Phil Renkert
    % Association: University of Illionis at Urbana-Champaign
    % Contact: aksland2@illinois.edu and renkert2@illinois.edu
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % potential improvements:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    properties
        Type extrinsicProp_Types = "Abstract" % Property type, i.e. mass, price, etc...
    end
    
    methods
        function obj = extrinsicProp(type, val, varargin)
            if isa(type, 'string') || isa(type, 'char')
                type = extrinsicProp_Types(type);
            elseif ~isa(type, 'extrinsicProp_Types')
                error('Type must be a valid extrinsic prop type');
            end
            
            obj = obj@compParam(string(type), val, 'AutoRename', true, varargin{:});
            obj.Type = type;
        end
    end
    
    methods (Sealed)
        function resProps = Combine(obj_array)
            % extrinsicProp.Combine combines like types in array of extrinsicProps 
            % into the resulting system extrinsicProps
            
            types = [obj_array.Type];
            unique_types = unique(types);
            
            resProps = extrinsicProp.empty();
            for i = 1:numel(unique_types) % For each unique type in obj_array
                prop_id = (unique_types(i) == types); % Get indices of all properties of that type
                breakpoints = obj_array(prop_id); % Breakpoints for Dependent compParam
                type = unique_types(i);
                units = [breakpoints.Unit];
                unit = units(1);
                if any(units ~= unit)
                    warning('Combining extrinsicProps with different units!')
                end
                F = @(varargin) type.combFunc([varargin{:}]);
                resProps(i,1) = extrinsicProp(unique_types(i), NaN,...
                    'Unit', unit,...
                    'Dependent',true,...
                    'DependentBreakpoints', breakpoints,...
                    'DependentFunction', F); % Assign aggregate prop value into new system extrinsicProp
            end     
        end
        
        function [val, i] = getProp(obj_array, type, parent)
            % Get value of extrinsicProp of Type type from object array of extrinsicProps
            i = [obj_array.Type] == type;
            if nargin == 3
                i = i & ([obj_array.Parent] == parent);
            end
            val = vertcat(obj_array(i));
        end
    end
end

