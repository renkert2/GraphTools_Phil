classdef Limits < handle
    % Limits is a class that provides bounding information in the Graph 
    % Model Toolbox. Limits will convert upper and lower bound box
    % constraint information into an H-representation. Instatiate this as
    % an empty object or by:
    % 
    % L = Limit(bnd1,bnd2) where bnd1 is one bound for the parent object
    % and bnd2 is a second bound for the parent object. Upper and lower
    % bounds are assigned automatically
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Contributors: Christopher T. Aksland and Phil Renkert
    % Association: University of Illionis at Urbana-Champaign
    % Contact: aksland2@illinois.edu and renkert2@illinois.edu
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % potential improvements:
    % integrate with the MPT toolbox
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    properties
        UpperLimit (1,1) double = inf 
        LowerLimit (1,1) double = -inf
    end
    
    properties (Dependent)
        A (:,:) % H-rep such that {x|Ax<=b}
        b (:,1) % H-rep such that {x|Ax<=b}
    end
    
    methods
        
        function obj = Limits(varargin) % constructor method
            if nargin == 0
                
            elseif nargin == 2
                obj.UpperLimit = max([varargin{:}]);
                obj.LowerLimit = min([varargin{:}]);
                if obj.UpperLimit == obj.LowerLimit % throw a warning if the limits are equal
                    warning('Object has the same upper and lower limit')
                end
            else % throw an error if too much or little information is provided at object creation
                error('Define both upper and lower limits at object creation')
            end
        end
        
        function A = get.A(obj) % pull the A matrix in the H-rep
            H = getHRep(obj);
            A = H(:,1);
            
        end
        
        function b = get.b(obj) % pull the b matrix in the H-rep
            H = getHRep(obj);
            b = H(:,2);
        end
        
        function H = getHRep(obj) % build the model H-representation
            H = zeros(2,2);
            H(:,1) = [1;-1];
            H(:,2) = [obj.UpperLimit;-obj.LowerLimit];
                      
        end
        
        
        
    end
    
end