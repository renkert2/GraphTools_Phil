classdef compParam < handle
    % compParams are used to define component parameters.  symParams can
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
        Sym string
        Value double
        Tunable logical = false
        AutoRename logical = false % Automatically append name of Parent element to end of Sym and Sym_
        Description string
        Unit
        
        Parent SystemElement
    end
    
    properties (SetAccess = private, Hidden = true)
        Sym_ sym
    end
    
    methods
        function obj = compParam(sym_arg, def, opts)
            arguments
                sym_arg string  = ""
                def double = 0
                opts.Tunable = false
                opts.Description = ""
                opts.Unit = ""
            end
            
            obj.Sym = sym_arg; 
            obj.Value = def;
            obj.Tunable = opts.Tunable;
            obj.Description = opts.Description;
            obj.Unit = opts.Unit;
        end
        
        function x = pop(obj)
            % pop() returns symbolic object or default value depending on 
            % tunable property.  Detaches value from compParam handle object
            if obj.Tunable
                x = obj.Sym_;
            else
                x = obj.Value;
            end
        end
        
        function d = double(obj)
            d = pop(obj);
        end
        
        function x = sym(obj)
            x = sym(pop(obj));
        end
        
        function set.Sym(obj, symarg)   
            if obj.AutoRename && ~isempty(obj.Parent)
                symarg = obj.appendParentName(symarg);
            end
            
            if obj.Sym ~= symarg
                obj.Sym = symarg;
                obj.Sym_ = symarg.empty();
            end
        end
        
        function x = get.Sym_(obj)
            % Get method keeps us from having to use the symbolic toolbox more than necessary
            if isempty(obj.Sym_)
                obj.setSym_();
            end
            x = obj.Sym_;
        end
        
        function set.Parent(obj, par)
            obj.Parent = par;
            obj.Sym = obj.Sym; % Rename sym with Parent name if AutoRename flag is true
        end
    end
    
    %% Private
    methods (Hidden = true)
        function setSym_(obj)
            obj.Sym_ = sym(obj.Sym, ["real", "positive"]); % Add real, positive assumptions to all sym_param symbolic variables
        end
        
        function str_out = appendParentName(obj, str_in)
            p_name = obj.Parent.Name;
            p_name = regexprep(p_name, '[\W]', '_');
            str_out = str_in+"__"+p_name;
        end
    end
    %% Overload!
    methods
        function x = times(varargin)
            x = compParam.processOverload(@times, varargin);
        end
        
        function x = mtimes(varargin)
            x = compParam.processOverload(@mtimes, varargin);
        end
        
        function x = plus(varargin)
            x = compParam.processOverload(@plus, varargin);
        end
        
        function x = uplus(varargin)
            x = compParam.processOverload(@uplus, varargin);
        end
        
        function x = minus(varargin)
            x = compParam.processOverload(@minus, varargin);
        end
        
        function x = rdivide(varargin)
            x = compParam.processOverload(@rdivide, varargin);
        end
        
        function x = mrdivide(varargin)
            x = compParam.processOverload(@mrdivide, varargin);
        end
    end
    methods (Static)
        function x = processOverload(func, args)
            for i = 1:numel(args)
                if isa(args{i}, 'compParam')
                    args{i} = pop(args{i});
                end
            end
            
            x = func(args{:});
        end            
    end
end

