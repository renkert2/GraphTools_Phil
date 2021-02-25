classdef Limits < handle
    %   Limits class
    %   This class is used to define limits for aspects of a graph.
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Author: Christopher T. Aksland 
    % Association: University of Illionis at Urbana-Champaign
    % Contact: aksland2@illinois.edu
    % Revision History:
    % 2/1/2021 - Class Creation
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % potential improvements:
    % Probably should shift this to rely on MPT polyhedron objects as
    %   constraints. However, our implementation is probably more
    %   computationally lean.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    properties
        UpperLimit (1,1) double = inf
        LowerLimit (1,1) double = -inf
    end
    
    properties (Dependent)
        A % H-rep such that {x|Ax<=b}
        b % H-rep such that {x|Ax<=b}
    end
    
    methods
        
        function obj = Limits(varargin)
            if nargin == 0
                
            elseif nargin == 2
                obj.UpperLimit = max([varargin{:}]);
                obj.LowerLimit = min([varargin{:}]);
                if obj.UpperLimit == obj.LowerLimit
                    warning('Object has the same upper and lower limit')
                end
            else
                error('Define both upper and lower limits at object creation')
            end
        end
        
        function A = get.A(obj)
            H = getHRep(obj);
            A = H(:,1);
            
        end
        
        function b = get.b(obj)
            H = getHRep(obj);
            b = H(:,2);
        end
        
        function H = getHRep(obj)
            H = zeros(2,2);
            H(:,1) = [1;-1];
            H(:,2) = [obj.UpperLimit;-obj.LowerLimit];
                      
        end
        
        
        
    end
    
    
    
end