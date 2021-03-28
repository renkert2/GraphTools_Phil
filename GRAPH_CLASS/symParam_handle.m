classdef symParam_handle < handle
    % symParams are used to define symbolic design variables 
    % for design optimization or parameter studies.  symParams can
    % replace numeric parameters in a component definition.  When using 
    % symParams in a GraphModel, capacitance and powerflow coefficients 
    % become symbolic.  Also, CalcF and CalcG require an additional 
    % argument for with numeric values for the symbolic parameters.
    % 
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Contributors: Christopher T. Aksland and Phil Renkert
    % Association: University of Illionis at Urbana-Champaign
    % Contact: aksland2@illinois.edu and renkert2@illinois.edu
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % potential improvements:  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    properties
        Sym sym
        Default_Value double
    end
    
    methods
        function obj = symParam_handle(sym_arg, def, varargin)  
            % symParam Constructor:
            % - sym_arg: char or string of the name of the symbolic variable, i.e. "R" or "J"
            % - def: Default value of the symParam, double
            % - varargin: passes additional options to the sym function
            
            if nargin <= 1
                def = 0;
            end
            
            if nargin == 0
                sym_arg = [];
            end
            
            obj.Sym = sym(sym_arg, ["real", "positive", varargin{:}]); % Add real, positive assumptions to all sym_param symbolic variables
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

