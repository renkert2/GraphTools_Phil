classdef compParam < handle & matlab.mixin.Heterogeneous
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
        Size (1,2) double = [1 1]
        Value double = NaN
        Tunable logical = false
        AutoRename logical = true % Automatically append name of Parent element to end of Sym and Sym_
        Description string = ""
        Unit string = ""
        Assumptions string = ["real", "positive"]
        
        Parent SystemElement
        
        % Functionality for Dependent Parameters.  Could also be made a subclass
        Dependent logical = false % Flag that states whether to reference call internal function or cached obj.Value
        DependentFunction % Function used to evaluate the dependent parameter value - must accept comma seperated list of inputs
        DependentBreakpoints (:,1) compParam % compParams on which obj depends.  These are the arguments for obj.DependenceFunction.  Leave empty if DependenceFunction takes no arguments
    end
    
    properties (SetAccess = private)
        SymID string 
        Sym_ sym
        DependentDefault logical
    end
    
    methods
        function obj = compParam(sym_arg, val, varargin)
            if nargin >= 1
                obj.Sym = sym_arg; 
            end
            if (nargin >= 2)
                obj.Size = size(val);
                obj.Value = val;
            end
            if nargin > 2
                my_inputparser(obj,varargin{:});
            end
            update(obj)
        end
        
        function setDependency(obj, fun, bkpnts)
            % Makes adding dependencies a little faster in construction
            % obj.Dependent = true;
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
        
        function set.Size(obj, size)
            obj.Size = size;
            obj.Sym_ = sym.empty();
        end
        
        function set.Value(obj, val)
            assert(all(size(val) == obj.Size), 'Size of Value must match obj.Size')
            obj.Value = val;
        end
        
        function set.Assumptions(obj, assumptions)
            obj.Assumptions = assumptions;
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
            obj.Sym = obj.Sym; % Call Set Sym set method to Rename sym with Parent name if AutoRename flag is true
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
                        if iscell(bkpntsPop)
                            x = obj.DependentFunction(bkpntsPop{:});
                        else
                            x = obj.DependentFunction(bkpntsPop);
                        end
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
    
        function i = find(obj, args1, args2)
            if nargin == 2
                % Single argument -> find by SymID
                i = arrayfun(@(sid) find(indexBy(obj, "SymID", sid)), string(args1));
            elseif nargin == 3
                % Two arguments -> find by Sym, Parent
                i = arrayfun(@(s,p) find(indexBy(obj, "Sym", s) & indexByParent(obj, "Name", p)), string(args1), string(args2));
            end
        end
        
        function o = get(obj, varargin)
            o = obj(find(obj, varargin{:}));
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
        
        function modified_objs = loadValues(obj_array, data)
            % Updates values to array of compParamValues or ComponentData
            % object
            % Returns an array of compParams whose Values were set
            if isa(data, 'ComponentData')
                data = vertcat(data.Data);
            end
            
            obj_array_size = [numel(obj_array), 1];
            cp_syms = vertcat(obj_array.Sym);
            cp_comps = reshape(parentTypes(obj_array), obj_array_size);
            cp_idcols = [cp_syms, cp_comps];
            
            cpv_syms = vertcat(data.Sym);
            cpv_comps = vertcat(data.Component);
            cpv_idcols = [cpv_syms, cpv_comps];
            
            i_combined_accum = false(obj_array_size); % Logical array tracks which objects matched a compParamValue
            for i = 1:numel(data)% Loop over each compParamValue
                cpv_id = cpv_idcols(i,:);
                
                % Ensure each cpv can be uniquely identified
                assert(sum(all(cpv_id == cpv_idcols,2)) == 1, "Each compValue entry in data must have a unique combination of Sym and Component");
                
                i_combined = all(cpv_id == cp_idcols,2);
                if any(i_combined)
                    assert(sum(i_combined) == 1, "Multiple compParams correspond to compParamValue %d", i);
                    obj = obj_array(i_combined);
                    
                    
                    unit_test = @(s) ~(isempty(s) || s == "");
                    cpv_unit = data(i).Unit;
                    obj_unit = obj.Unit;
                    if unit_test(cpv_unit) && unit_test(obj_unit)
                        assert(strcmpi(cpv_unit, obj_unit), "Incompatible Units in %s parameter %s: compParam has unit %s and compParamValue has unit %s",...
                            data(i).Component, data(i).Sym, obj_unit, cpv_unit);
                    end
                    
                    obj.Value = data(i).Value;
                    
                    i_combined_accum = i_combined_accum | i_combined;
                end
            end
            modified_objs = obj_array(i_combined_accum);
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
              
        function [f, f_mtlb] = matlabFunction(obj_array, sym, args, opts)
            % matlabFunction generates a matlabFunction from a syms 
            % object using the SymParams.  If sym is not symbolic, generate
            % a standard function that returns the constant value regardless
            % of the input
            % Args can be used to specify additional input arguments for
            % other symbolic variables not in the SymParams.  Elements of
            % cell array args are added as arguments after the first
            % argument for the sym param values
            
            % f automatically evaluates params.Value and substitutes them into the 
            % function so that they are hidden away.  Probably a better way to do this
            % but does work well.  
            arguments
                obj_array
                sym
                args cell = {}
                opts.FilterTunable logical = true
            end
            
            if isa(sym, 'sym')
                if opts.FilterTunable
                    params = obj_array(isTunable(obj_array));
                else
                    params = obj_array;
                end
                
                if ~isempty(params)
                    f_mtlb = matlabFunction(sym, 'Vars', [{params.Sym_}, args]); % params.Sym_ is a comma separated list
                    f = @(varargin) f_mtlb(params.Value, varargin{:}); % params.Value is a comma separated list
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
                opts.InlineArg = "$$"
            end
            
            N = numel(obj_array);
            s = string.empty(N,0);
            for i = 1:N
                obj = obj_array(i);
                s_temp = opts.InlineArg+obj.Sym+opts.InlineArg;
                if ~isempty(obj.Unit) && obj.Unit ~= "" && opts.UnitFlag
                    s_temp = s_temp + " "+"("+obj.Unit+")";
                end
                s(i,1) = s_temp;
            end
        end
        
        function storeDependentDefault(obj_array)
            for i = 1:numel(obj_array)
                obj = obj_array(i);
                obj.DependentDefault = obj.Dependent;
            end
        end
        
        function dependent_old = setDependent(obj_array, dependent_new)
            dependent_old = vertcat(obj_array.Dependent);
            for i = 1:numel(obj_array)
                obj = obj_array(i);
                if isscalar(dependent_new)
                    obj.Dependent = dependent_new;
                else
                    obj.Dependent = dependent_new(i);
                end
            end
        end
        
        function restoreDependentDefault(obj_array)
            for i = 1:numel(obj_array)
                obj = obj_array(i);
                if ~isempty(obj.DependentDefault)
                    obj.Dependent = obj.DependentDefault;
                end
            end
        end
        
        function tbl = table(obj_array, varargin)
            tbl = dispTable(obj_array, varargin{:});
        end
        
        function tbl = dispTable(obj_array, fields, opts)
            arguments
                obj_array
                fields string = ["Sym", "Value", "Unit", "Description", "Parent", "Dependent"]
                opts.LatexSym logical = false
                opts.LatexOpts cell = {}
            end
            field_vals_cell = cell(size(fields));
            for i = 1:numel(fields)
                field = fields(i);
                switch field
                    case "Value"
                        try
                            vals = vertcat(obj_array.Value);
                        catch
                            vals = strings(numel(obj_array),1);
                            for j = 1:numel(vals)
                                val = obj_array(j).Value;
                                if isnumeric(val) && isscalar(val)
                                    vals(j,1) = num2str(val);
                                else
                                    sz = size(val);
                                    type = class(val);
                                    vals(j,1) = sprintf("%dx%d %s", sz(1), sz(2), type);
                                end
                            end
                        end
                        field_vals = vals;
                    case "Sym"
                        if opts.LatexSym
                            sym = latex(obj_array, opts.LatexOpts{:});
                        else
                            sym = vertcat(obj_array.Sym);
                        end
                        field_vals = sym;
                    case "Parent"
                        field_vals = vertcat(vertcat(obj_array.Parent).Name);
                    otherwise
                        field_vals = vertcat(obj_array.(field));
                end
                field_vals_cell{i} = field_vals;
            end
            
            tbl = table(field_vals_cell{:},'VariableNames', fields);
            if nargout == 0
                disp(tbl);
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
    
    %% Static
    methods (Static)
        function cp = gatherObjectParams(obj)
            % Collects compParams assigned to object's properties into a single compParam array
            obj_meta = metaclass(obj);
            parent_obj = metaclass(SystemElement);
            parent_props = parent_obj.PropertyList;
            meta_props = setdiff(obj_meta.PropertyList, parent_props);
            
            % Make property filter
            meta_filter = ~[meta_props.Dependent]; % Initially exclude dependent props
            for i = 1:numel(meta_filter)
                if meta_filter(i)
                    meta_prop = meta_props(i);
                    v = meta_prop.Validation;
                    if isempty(v) || isempty(v.Class)
                        meta_filter(i) = false;
                    else
                        meta_cp = ?compParam;
                        v_class_flag = (v.Class == meta_cp) || ismember(meta_cp, v.Class.SuperclassList);
                        if ~v_class_flag
                            meta_filter(i) = false;
                        end
                    end
                end
            end
            
            meta_props = meta_props(meta_filter);
            
            cp = compParam.empty();
            for i = 1:numel(meta_props)
                prop = obj.(meta_props(i).Name);
                if isa(prop, 'compParam') && isscalar(prop)
                    cp(end+1,1) = prop;
                end
            end 
        end
    end
    
    %% Private
    methods (Hidden = true)
        function setSym_(obj)
            if all(obj.Size == [1 1])
                obj.Sym_ = sym(obj.SymID, obj.Assumptions);
            else
                obj.Sym_ = sym(obj.SymID, obj.Size, obj.Assumptions);
            end
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

