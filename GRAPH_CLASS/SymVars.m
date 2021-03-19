classdef SymVars
    % SymVars stores the symbolic variables x,u,d for the GeneralModel Class.  Replaces SymVars struct() in Model
    properties
        x (:,1) sym
        d (:,1) sym
        u (:,1) sym
    end
    
    properties (Hidden)
        x_full (:,1) sym % Used by GraphModel
    end
    
    methods
        function obj = SymVars(opts)
            arguments
                opts.Nx double = []
                opts.Nd double = []
                opts.Nu double = []
            end
            
            if ~isempty(opts.Nx)
                x = obj.genSymVars('x%d', opts.Nx); % Symbolic Vars for Dynamic States
                obj.x = x;
            end
            if ~isempty(opts.Nd)
                d = obj.genSymVars('d%d', opts.Nd); % Symbolic Vars for Disturbances
                obj.d = d;
            end
            if ~isempty(opts.Nu)
                u = obj.genSymVars('u%d', opts.Nu); % Symbolic Vars for Inputs
                obj.u = u;
            end
        end
    end
    
    methods (Static)
        function s = genSymVars(patt, n)
            if n == 1
                s = sym(sprintf(patt,1));
            elseif n>1
                s = sym(patt, [n 1]);
            else
                s = sym.empty();
            end
        end
    end
end