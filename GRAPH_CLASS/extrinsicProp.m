classdef extrinsicProp
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
        Value {mustBeA(Value, ["double", "sym"])} % Numeric or Symbolic value associated with the property
    end
    
    methods
        function obj = extrinsicProp(type, val)
            if nargin == 2
                obj.Type = type;
                obj.Value = val;
            end
        end
        
        function resProps = Combine(obj_array)
            % extrinsicProp.Combine combines like types in array of extrinsicProps 
            % into the resulting system extrinsicProps
            
            types = [obj_array.Type];
            unique_types = unique(types);
            
            resProps = extrinsicProp.empty();
            for i = 1:numel(unique_types) % For each unique type in obj_array
                prop_id = (unique_types(i) == types); % Get indices of all properties of that type
                val = unique_types(i).combFunc([obj_array(prop_id).Value]); % Call the combination function of the type on the values of the property of that type
                resProps(i) = extrinsicProp(unique_types(i), val); % Assign aggregate prop value into new system extrinsicProp
            end     
        end
        
        function [val, i] = getProp(obj_array, type)
            % Get value of extrinsicProp of Type type from object array of extrinsicProps
            i = [obj_array.Type] == type;
            val = [obj_array(i).Value];
        end
    end
end

