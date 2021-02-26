classdef symParam < sym
    % symParam ...
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
        Default_Value double
    end
    
    methods
        function obj = symParam(sym_arg, def, varargin)            
            if nargin <= 1
                def = 0;
            end
            
            if nargin == 0
                sym_arg = [];
            end
            
            obj@sym(sym_arg, ["real", "positive", varargin{:}]);
            obj.Default_Value = def;

        end
        
        function d = double(obj)
            d = obj.Default_Value;
        end
        
        function d = getDefault(obj)
            d = obj.Default_Value;
        end
    end
end

