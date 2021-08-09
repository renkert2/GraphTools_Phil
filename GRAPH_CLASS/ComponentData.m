classdef ComponentData
    %Value class encapsulating data corresponding to a physical component
    % It contains an array of compParamValues in its Data property that are
    % used to replace the values of system element parameters.  
    
    properties
        Component string % Corresponds to Component class name
        Make string
        Model string
        SKU string
        Description string
        
        Data (:,1) compParamValue
    end
    
    methods
        function obj = ComponentData(comp, make, model, data)
            if nargin > 0
                obj.Component = comp;
                obj.Make = make;
                obj.Model = model;
                obj.Data = data(:);

                for i = 1:numel(obj.Data)
                    obj.Data(i).Component = obj.Component;
                end
            end
        end
        
        function [CD, unique_comps, N_comps] = separateComponents(obj_array, des_vars)
            arguments
                obj_array
                des_vars compParamValue = compParamValue.empty() % Optional, only includes components in des_vars
            end
            
            if nargin >1
                unique_comps = intersect(unique([des_vars.Component],'stable'), unique([obj_array.Component],'stable'), 'stable');
            else
                unique_comps = unique([obj_array.Component]);
            end
            
            N_unique_comps = numel(unique_comps);
            N_comps = zeros(1,N_unique_comps);
            CD = cell.empty(0,N_unique_comps); % Cell Array containing component data filtered by component
            for i = 1:N_unique_comps
                comp = unique_comps(i);
                [cd_comp, ~] = filterComponent(obj_array, comp);
                CD{1,i} = cd_comp;
                N_comps(i) = numel(cd_comp);
            end
        end
        
        function [filtered_cd] = filterBounds(obj_array, des_vars, LB, UB)
            % Regect entries in the catalog that violate bounds
            arguments
                obj_array
                des_vars compParamValue % Just used to find corresponding parameters, values not used
                LB (:,1) double
                UB (:,1) double
            end
            assert(numel(LB) == numel(des_vars), "LB must be same size as des_vars");
            assert(numel(UB) == numel(des_vars), "UB must be same size as des_vars")
            
            I = false(size(obj_array));
            for i = 1:numel(obj_array)
                candidate = obj_array(i);
                candidate_cpv = candidate.Data;
                [cpv_pairs, I_pairs] = getPairs(des_vars, candidate_cpv);
                I(i) = isInBounds(cpv_pairs(:,2), LB(I_pairs(:,1)), UB(I_pairs(:,1)));
            end
            filtered_cd = obj_array(I);
        end
        
        function [CD_sorted, D_sorted, unique_comps] = filterNearest(obj_array, target, N_max, opts)
            arguments 
                obj_array
                target compParamValue
                N_max = inf
                opts.DistanceMode string = "Norm"
                opts.Gradient double = []
                opts.Hessian double = []
                opts.Weights = []
                opts.Tolerance double = []
            end
            
            unique_comps = intersect(unique([target.Component],'stable'), unique([obj_array.Component],'stable'), 'stable');
            N_unique_comps = numel(unique_comps);
            
            CD_sorted = cell.empty(0,N_unique_comps); % Cell Array containing component data filtered by component
            D_sorted = cell.empty(0,N_unique_comps); % Cell array containing distances for each component
            
            for i = 1:N_unique_comps
                comp = unique_comps(i);
                
                if isa(N_max, 'struct')
                    n_max = N_max.(comp);
                elseif isscalar(N_max)
                    n_max = N_max;
                else
                    error("N_max argument must be scalar double or a struct where each field is a component with corresponding value N_max")
                end
                
                [cd_comp, ~] = filterComponent(obj_array, comp);
                [target_comp, target_I] = filterComponent(target, comp);

                weight = [];
                grad = [];
                hessian = [];
                tolerance = [];
                map = [];
                if any(target_I)
                    switch opts.DistanceMode
                        case "Norm"
                            % do nothing
                        case "WeightedNorm"
                            if isa(opts.Weights, 'struct')
                                weight = opts.Weights.(comp);
                            else
                                weight = opts.Weights(target_I);
                            end
                        case "TaylorSeries"
                            if ~isempty(opts.Gradient)
                                grad = opts.Gradient;
                            end
                            if ~isempty(opts.Hessian)
                                hessian = opts.Hessian;
                            end
                            map = find(target_I);
                            if ~isempty(opts.Tolerance)
                                if isscalar(opts.Tolerance)
                                    tolerance = opts.Tolerance;
                                else
                                    tolerance = opts.Tolerance(target_I);
                                end
                            end
                    end
                end

                [cd_comp_sorted, d_sorted] = processComp(target_comp, cd_comp, n_max, weight, grad, hessian, map, tolerance);
                
                CD_sorted{1,i} = cd_comp_sorted;
                D_sorted{1,i} = d_sorted;
            end
            
            if N_unique_comps == 1
                CD_sorted = CD_sorted{:};
                D_sorted = D_sorted{:};
            end
            
            function [cd_sorted,d_sorted] = processComp(target, cd, n_max, weights, grad, hessian, map, tolerance)
                N_comp = numel(cd);
                d = zeros(N_comp,1);
                for j = 1:N_comp
                    d(j) = distance(target, cd(j).Data, 'Mode', opts.DistanceMode,'Weights', weights, 'Gradient', grad, 'Hessian', hessian, 'Map', map, 'Tolerance', tolerance);
                end
                d = d(~isnan(d));
                [~, i_sorted] = sort(d, 'ascend');
                i_sorted = i_sorted(1:min(numel(i_sorted), n_max)); % Extract no more than N_max values
                
                cd_sorted = cd(i_sorted);
                d_sorted = d(i_sorted);
            end
        end
        
        function [cd_filtered, I] = filterComponent(obj_array, component)
           I = vertcat(obj_array.Component) == component;
           cd_filtered = obj_array(I); 
        end
        
        function [combinations, I] = combinations(varargin)
            % Enumerates all component combinations that can be formed by 
            % selecting one component out of each ComponentData array
            % argument.  The output is a PxN matrix, where P is the product
            % of the lengths of the arguments and N is the number of
            % arguments.  See allcomb documentation on file exchange for
            % more details
            arg_chk_1 = cellfun(@(x) isa(x, 'ComponentData'), varargin);
            arg_chk_2 = cellfun(@(x) min(size(x)) == 1, varargin);
            assert(all(arg_chk_1), 'All input arguments must be ComponentData vectors');
            assert(all(arg_chk_2), "All arguments must have at least one element");
            
            lens = cellfun(@(v) 1:numel(v), varargin, 'UniformOutput', false);
            I = allcomb(lens{:});
            combinations = allcomb(varargin{:});
        end
        
        function [tbl,param_table] = table(obj_array)
            comps = unique([obj_array.Component]);
            N = numel(comps);
            if N == 1
                [tbl, param_table] = makeTable(obj_array);
            else
                tbls = cell.empty(0,N);
                fields = cell.empty(0,N);
                param_table = cell.empty(0,N);
                common_fields = string.empty();
                for c = 1:N
                    [tbls{c}, param_table{c}] = makeTable(filterComponent(obj_array, comps(c)));
                    fields{c} = string(tbls{c}.Properties.VariableNames);
                    if c == 1
                        common_fields = fields{c};
                    else
                        common_fields = intersect(common_fields, fields{c});
                    end
                end
                common_tbls = cell.empty(0,N);
                for c = 1:N
                    FI = arrayfun(@(f) ismember(f,common_fields), fields{c});
                    common_tbls{c} = tbls{c}(:,FI);
                end
                tbl = [tbls, {vertcat(common_tbls{:})}];
            end

            function [tbl, param_table] = makeTable(obj_array)
                pdat_all = vertcat(obj_array.Data);
                [params, I] = unique(vertcat(pdat_all.Sym));
                unique_pdat = pdat_all(I);
                
                param_table = table(unique_pdat);
                pvars = string(param_table.Properties.VariableNames);
                pvi = pvars ~= "Value";
                param_table = param_table(:,pvi);
                
                comp_fields = ["Make", "Model", "SKU"];
                varnames = [comp_fields, params'];
                vartypes = [repmat("string",1,numel(comp_fields)) repmat("double",1,numel(params))];
                tbl = table('Size', [numel(obj_array), numel(varnames)], 'VariableTypes', vartypes, 'VariableNames', varnames);
                for i = 1:numel(obj_array)
                    for f = comp_fields
                        val = obj_array(i).(f);
                        if isempty(val)
                            val = "";
                        end
                        tbl(i,:).(f) = val;
                    end
                    dat = obj_array(i).Data;
                    for j = 1:numel(params)
                        f = params(j);
                        cpv = filterSym(dat, f);
                        if isempty(cpv)
                            val = NaN;
                        else
                            N_cpv = numel(cpv);
                            if N_cpv ~= 1
                                error("Multiple compParamVals of same name")
                            else
                                val = cpv.Value;
                            end
                        end
                        tbl(i,:).(f) = val;
                    end
                end
            end
        end
        
        function dispTable(obj_array, opts)
            arguments
               obj_array
               opts.ParamTable logical = false
            end
            
            comps = unique([obj_array.Component]);
            for c = comps
                fprintf("Component: %s \n", c)
                [t,tp] = table(filterComponent(obj_array, c));
                disp(t);
                if opts.ParamTable
                    disp(tp);
                end
                disp(newline);
            end
        end
        
        function tbl = summaryTable(obj_array)
            comp_fields = ["Component", "Make", "Model", "SKU", "Description"];
            
            tbldat = cell.empty(0,numel(comp_fields));
            for i = 1:numel(comp_fields)
                S = strings(numel(obj_array),1);
                for j = 1:numel(obj_array)
                    s = obj_array(j).(comp_fields(i));
                    if ~isempty(s)
                        S(j) = s;
                    end
                end
                tbldat{1,i} = S;
            end

            tbl = table(tbldat{:}, 'VariableNames', comp_fields);
        end
    end
    
    methods (Static)
        function obj_array = importFromFile(file)
            arguments
                file string
            end
            
            data = readcell(file, 'TextType', 'string', 'MissingRule', 'fill');
            
            % Breakpoints: Component ; empty ; Parameter Data; empty ; Component Data
            first_col = data(:,1);
            [ir_comp,ic_comp] = find(cellfun(@(x) matchText(x,"Component"), data));
            [ir_param, ~] = find(cellfun(@(x) matchText(x,"Parameter"), data));
            [ir_data, ic_data] = find(cellfun(@(x) matchText(x,"Data"), data));
            
            % Read Component
            component = data{ir_comp,ic_comp+1};
            
            % Read Parameters
            param_props = [first_col{ir_param:ir_data-1}];
            param_data = data(ir_param:ir_data-1,2:end);
            missing_cols = all(cellfun(@ismissing, param_data),1);
            param_data = param_data(:,~missing_cols)';
            param_table = cell2table(param_data, 'VariableNames', param_props);
            params = param_table{:,1}';
            
            % Component Properties
            comp_props = rmmissing([data{ir_data+1,ic_data:end}]);
            
            % Get Data
            varnames = [comp_props params];
            data_table = cell2table(data(ir_data+2:end,:),'VariableNames',varnames);
            
            % Make Template for compParamValue array
            data_template = compParamValue.importFromTable(param_table);
            
            % Create ComponentData Array
            N = height(data_table);
            obj_array = ComponentData.empty(N,0);
            
            for i = 1:N
                cd = ComponentData();
                cd.Component = component;
                for j = 1:numel(comp_props)
                    field = comp_props(j);
                    cd.(field) = data_table(i,:).(field);
                end
                
                cpvdata = data_template;
                for j = 1:numel(params)
                    cpvdata(j).Value = data_table(i,:).(params(j));
                    cpvdata(j).Component = component;
                end
                cd.Data = cpvdata;
                obj_array(i,1) = cd;
            end
            
            % Helper Functions
            function l = matchText(str,patt)
                if isa(str, 'string') || isa(str, 'char')
                    l = contains(str,patt, 'IgnoreCase', true);
                else
                    l = false;
                end
            end
        end
        function obj_array = importFromJSON(file)
            arguments
                file string
            end
            raw = fileread(file);
            json = jsondecode(raw); % Struct array, each element is a component
            obj_array = ComponentData.importFromStruct(json);
        end
        function obj_array = importFromStruct(s)
            
            N = numel(s);
            obj_array = ComponentData.empty(N,0);
            
            for i = 1:N
                if isa(s,'cell')
                    comp_struct = s{i};
                elseif isa(s,'struct')
                    comp_struct = s(i);
                else
                    error('Input argument must be struct array or cell array of structs')
                end
                dat = comp_struct.Data;
                cpv = compParamValue.importFromStruct(dat);
                for j = 1:numel(cpv)
                    cpv(j).Component = comp_struct.Component;
                end
                comp_struct.Data = cpv;
                
                cd = ComponentData();
                field_names = string(fields(comp_struct))';
                for f = field_names
                    cd.(f) = comp_struct.(f);
                end
                    
                obj_array(i,1) = cd;
            end   
        end
    end
end

