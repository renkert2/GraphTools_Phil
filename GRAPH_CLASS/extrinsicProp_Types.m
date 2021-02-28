classdef extrinsicProp_Types
    % extrinsicProp_Types is an enumeration class that stores information 
    % pertaining to the different types of extrinsicProps, i.e. mass and 
    % price.  Each type is associated with a combination function (usually
    % @sum), and a numeric id that allows us to define a sort method for the 
    % class.  
    % A class must have a sort method to call unique() on the class
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Contributors: Christopher T. Aksland and Phil Renkert
    % Association: University of Illionis at Urbana-Champaign
    % Contact: aksland2@illinois.edu and renkert2@illinois.edu
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % potential improvements:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    enumeration
        Abstract (@sum, 1)
        Mass (@sum, 2)
        Cost (@sum, 3)
    end
    
    properties
        combFunc function_handle % function used to map enumeration values from Component -> System
        id uint8 % Used for functions like 'Sort'; Required to use @unique
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
