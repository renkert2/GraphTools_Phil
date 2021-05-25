classdef compParam < handle & matlab.mixin.Heterogeneous & matlab.mixin.CustomDisplay
    % compParams are used to define component parameters.  symParams can
    % replace numeric parameters in a component definition.
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
        Value double = NaN;
        Tunable logical = false
        AutoRename logical = true % Automatically append name of Parent element to end of Sym and Sym_
        Description string = ""
        Unit string = ""
        
        Parent SystemElement
        
        % Functionality for Dependent Parameters.  Could also be made a subclass
        Dependent logical = false % Flag that states whether to reference call internal function or cached obj.Value
        DependentFunction % Function used to evaluate the dependent parameter value - must accept comma seperated list of inputs
        DependentBreakpoints (:,1) compParam % compParams on which obj depends.  These are the arguments for obj.DependenceFunction.  Leave empty if DependenceFunction takes no arguments
    end
    
    properties (SetAccess = private, Hidden = true)
        SymID string 
        Sym_ sym
    end
    
    methods
        function obj = compParam(sym_arg, val, varargin)
            if nargin >= 1
                obj.Sym = sym_arg; 
            end
            if (nargin >= 2)
                obj.Value = val;
            end
            if nargin > 2
                my_inputparser(obj,varargin{:});
            end
            update(obj)
        end
        
        function setDependency(obj, fun, bkpnts)
            % Makes adding dependencies a little faster in construction
            obj.Dependent = true;
            obj.DependentFunction = fun;
            obj.DependentBreakpoints = bkpnts;
            obj.update();
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
                x = popScalar(obj_array);
            else
                x = cell.empty(numel(obj_array),0);
                for i = 1:numel(obj_array)
                    obj = obj_array(i);
                    x{i} = popScalar(obj);
                end
            end
            
            function x = popScalar(obj)
                if obj.Tunable
                    x = obj.Sym_;
                else
                    if ~obj.Dependent
                        x = obj.Value;
                    else
                        bkpntsPop = pop(obj.DependentBreakpoints); % Cell array of syms and/or doubles
                        x = obj.DependentFunction(bkpntsPop{:});
                    end
                end
            end
        end
        
        function update(obj_array)
            % Use update() method to update the value field of all dependent compParams
            % in a compParam object array.  This will recursively update nested dependent
            % parameters by updating all of the compParams in DependentBreakpoints.
            for i = 1:numel(obj_array)
                obj = obj_array(i);
                if obj.Dependent
                    obj.DependentBreakpoints.update();
                    calcDependentValue(obj);
                end
            end
        end
    
        function i = find(obj, args)
            i = arrayfun(@(x) find(indexBy(obj, "SymID", x)), string(args));
        end
        
        function o = get(obj, args)
            o = obj(find(obj, args));
        end
        
        function cpvals = getValues(obj_array)
            % 'Freezes' value of compParams and exports a compParamValue
            % array that can be reloaded later.  
            for i = 1:numel(obj_array)
                param = obj_array(i);
                pval = compParamValue(param.Sym, 'Value', param.Value);
                pval.Unit = param.Unit;
                pval.Description = param.Description;
                if ~isempty(param.Parent)
                    pval.Component = class(param.Parent);
                end
                cpvals(i,1) = pval;
            end
        end
        
        function loadValues(obj_array, data)
            % Updates values to array of compParamValues or ComponentData
            % object
            if isa(data, 'ComponentData')
                data = data.Data;
            end
            
            comps = parentTypes(obj_array);
            syms = reshape(vertcat(obj_array.Sym), size(obj_array));
            
            for i = 1:numel(data)
                i_comp = data(i).Component == comps;
                i_sym = data(i).Sym == syms;
                i_combined = i_comp & i_sym;
                if any(i_combined)
                    objs = obj_array(i_combined);
                    
                    unit = data(i).Unit;
                    if ~isempty(unit) || unit ~= ""
                        objs_units = vertcat(objs.Unit);
                        assert(all(objs_units == unit), "Incompatible Units in %s parameter %s", data(i).Component, data(i).Sym);
                    end
                    
                    for j = 1:numel(objs)
                        objs(j).Value = data(i).Value;
                    end
                end
            end
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
        
        function s = latex(obj_array, opts)
            arguments
                obj_array
                opts.UnitFlag = true
            end
            
            N = numel(obj_array);
            s = string.empty(N,0);
            for i = 1:N
                obj = obj_array(i);
                s_temp = "$$"+obj.Sym+"$$";
                if ~isempty(obj.Unit) && obj.Unit ~= "" && opts.UnitFlag
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
        
        function s = parentTypes(obj_array) 
            s = string.empty(numel(obj_array), 0);
            for i = 1:numel(obj_array)
                s(i) = class(obj_array(i).Parent);
            end
            s = reshape(s,size(obj_array));
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
        
        function calcDependentValue(obj)
            if ~isempty(obj.DependentFunction)
                obj.Value = obj.DependentFunction(obj.DependentBreakpoints.Value);
            end
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

