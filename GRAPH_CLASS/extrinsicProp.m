classdef extrinsicProp
    % extrinsicProp ...
    % @Phil
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Contributors: Christopher T. Aksland and Phil Renkert
    % Association: University of Illionis at Urbana-Champaign
    % Contact: aksland2@illinois.edu and renkert2@illinois.edu
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % potential improvements:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    properties
        Type extrinsicProp_Types = "Abstract"
        Value {mustBeA(Value, ["double", "sym"])}
    end
    
    methods
        function obj = extrinsicProp(type, val)
            if nargin == 2
                obj.Type = type;
                obj.Value = val;
            end
        end
        
        function resProps = Combine(obj_array)
            types = [obj_array.Type];
            unique_types = unique(types);
            
            resProps = extrinsicProp.empty();
            for i = 1:numel(unique_types)
                prop_id = (unique_types(i) == types);
                val = unique_types(i).combFunc([obj_array(prop_id).Value]);
                resProps(i) = extrinsicProp(unique_types(i), val);
            end     
        end
    end
end

