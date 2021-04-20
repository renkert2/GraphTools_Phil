classdef compParam < handle & matlab.mixin.Heterogeneous & matlab.mixin.CustomDisplay
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
        Value
        Tunable logical = false
        AutoRename logical = false % Automatically append name of Parent element to end of Sym and Sym_
        Description string
        Unit
        
        Parent SystemElement
    end
    
    properties (SetAccess = private, Hidden = true)
        SymID string 
        Sym_ sym
    end
    
    methods
        function obj = compParam(sym_arg, val, opts)
            arguments
                sym_arg string
                val
                opts.Tunable logical = false
                opts.AutoRename logical = false
                opts.Description = ""
                opts.Unit = ""
            end
            
            obj.Sym = sym_arg; 
            obj.Value = val;
            obj.Tunable = opts.Tunable;
            obj.AutoRename = opts.AutoRename;
            obj.Description = opts.Description;
            obj.Unit = opts.Unit;
        end
        
        function d = double(obj)
            d = pop(obj);
        end
        
        function x = sym(obj)
            x = sym(pop(obj));
        end
        
        function set.Sym(obj, symarg)  
            obj.Sym = symarg;
            if obj.AutoRename && ~isempty(obj.Parent)
                symsubs = obj.appendParentName(symarg);
                obj.SymID = symsubs;
            else
                obj.SymID = obj.Sym;
            end
            
            obj.Sym_ = sym.empty();
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
    
    %% Sealed Methods
    % Methods for heterogeneous object arrays of compParams
    
    methods (Sealed)  
        function x = pop(obj_array)
            % pop() returns symbolic object or default value depending on
            % tunable property.  Detaches value from compParam handle object
            if isscalar(obj_array)
                if obj_array.Tunable
                    x = obj_array.Sym_;
                else
                    x = obj_array.Value;
                end
            else
                x = sym.empty(numel(obj_array),0);
                for i = 1:numel(obj_array)
                    obj = obj_array(i);
                    if obj.Tunable
                        x(i) = obj.Sym_;
                    else
                        x(i) = obj.Value;
                    end
                end
            end
        end
    
        function i = find(obj, args)
            i = arrayfun(@(x) find(indexBy(obj, "SymID", x)), string(args));
        end
        
        function o = get(obj, args)
            o = obj(find(obj, args));
        end
        
        function exps = extrinsicProps(obj_array)
            i = isaArrayFun(obj_array, 'extrinsicProp');
            exps = obj_array(i);
        end
        
        function i = isTunable(obj_array)
            i = vertcat(obj_array.Tunable);
        end
        
        function [syms, tunable] = tunableSyms(obj_array)
            tunable = obj_array(isTunable(obj_array));
            syms = vertcat(tunable.Sym_);
        end
        
        function [vals, tunable] = tunableVals(obj_array)
            tunable = obj_array(isTunable(obj_array));
            vals = vertcat(tunable.Value);
        end
        
        function cp = unique(obj_array)
            syms = vertcat(obj_array.SymID);
            [~, ia, ~] = unique(syms);
            cp = obj_array(ia);
        end
              
        function [f, f_mtlb] = matlabFunction(obj_array, sym, args)
            % matlabFunction generates a matlabFunction from a syms 
            % object using the SymParams.  If sym is not symbolic, generate
            % a standard function that returns the constant value regardless
            % of the input
            % Args can be used to specify additional input arguments for
            % other symbolic variables not in the SymParams.  Elements of
            % cell array args are added as arguments after the first
            % argument for the sym param values
            arguments
                obj_array
                sym
                args cell = {}
            end
            
            if isa(sym, 'sym')
                tunables = tunableSyms(obj_array);
                if ~isempty(tunables)
                    f_mtlb = matlabFunction(sym, 'Vars', [{tunables}, args]);
                    f = @(varargin) f_mtlb(tunableVals(obj_array), varargin{:});
                else
                    f_mtlb = matlabFunction(sym, 'Vars', args);
                    f = f_mtlb;
                end
                    
            else
                f_mtlb = @(varargin) sym;
                f = f_mtlb;
            end
        end
        
        function s = exportStruct(obj_array, opts)
            arguments
               obj_array
               opts.ExportUnits = false
            end
            
            s = struct();
            for i = 1:numel(obj_array)
                obj = obj_array(i);
                val = obj.Value;
                if isa(val, 'sym')
                    val = double(subs(val, obj_array.tunableSyms, obj_array.tunableVals));
                end
                
                if opts.ExportUnits
                   struct_val = struct();
                   struct_val.Value = val;
                   struct_val.Unit = obj.Unit;
                else
                    struct_val = val;
                end
                
                comp_name = obj.Parent.Name;
                sym_name = genvarname(obj.Sym);
                if ~isempty(comp_name)
                    comp_name = genvarname(comp_name);
                    s.(comp_name).(sym_name) = struct_val;
                else
                    s.(sym_name) = struct_val;
                end    
            end
        end
        
        function s = latex(obj_array)
            N = numel(obj_array);
            s = string.empty(N,0);
            for i = 1:N
                obj = obj_array(i);
                s_temp = "$$"+obj.Sym+"$$";
                if ~isempty(obj.Unit) && obj.Unit ~= ""
                    s_temp = s_temp + " "+"("+obj.Unit+")";
                end
                s(i,1) = s_temp;
            end
        end
    end
    
    methods (Sealed, Hidden = true)
        function i = indexBy(obj_array, props, vals)
            if isscalar(props)
                i = vertcat(obj_array.(props)) == vals;
            else
                i_array = zeros(numel(obj_array), numel(props));
                for j = 1:numel(props)
                    i_array(:,j) = vertcat(obj_array.(props(j))) == vals(j);
                end
                i = all(i_array,2);
            end
        end
        
        function i = indexByParent(obj_array, props, vals)
            parents = vertcat(obj_array.Parent);
            if isscalar(props)
                i = vertcat(parents.(props)) == vals;
            else
                i_array = zeros(numel(parents), numel(props));
                for j = 1:numel(props)
                    i_array(:,j) = vertcat(parents.(props(j))) == vals(j);
                end
                i = all(i_array,2);
            end
        end
        
        function obj_filt = filterBy(obj_array, props, vals)
            obj_filt = obj_array(indexBy(obj_array, props, vals));
        end
        
        function obj_filt = filterByParent(obj_array, props, vals)
            obj_filt = obj_array(indexByParent(obj_array, props, vals));
        end
        
        function X = filteredProps(obj, prop, filter)
            vec = vertcat(obj.(prop));
            X = vec(filter);
        end
        
        function i = isaArrayFun(obj_array, classname)
            i = arrayfun(@(x) builtin('isa',x,classname), obj_array);
        end
    end
    
        
    methods (Sealed, Access = protected)
        function displayNonScalarObject(obj_array)
            % Modifies method from Mixin.Custom Display
            tbl = table(vertcat(obj_array.Sym), string(vertcat(obj_array.Value)), vertcat(obj_array.Unit), vertcat(obj_array.Description), vertcat(vertcat(obj_array.Parent).Name),...
                'VariableNames', ["Sym", "Value", "Unit", "Description", "Parent"]);
            disp(tbl);
        end
    end
    
    %% Private
    methods (Hidden = true)
        function setSym_(obj)
            obj.Sym_ = sym(obj.SymID, ["real", "positive"]); % Add real, positive assumptions to all sym_param symbolic variables
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
        
        function x = power(varargin)
            x = compParam.processOverload(@power, varargin);
        end
        
        function x = mpower(varargin)
            x = compParam.processOverload(@mpower, varargin);
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

