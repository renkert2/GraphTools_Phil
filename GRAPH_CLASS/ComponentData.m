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
        
        function [CD_sorted, D_sorted, unique_comps] = filterNearest(obj_array, target, N_max, opts)
            arguments 
                obj_array
                target compParamValue
                N_max = inf
                opts.LB (:,1) double
                opts.UB (:,1) double
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
                
                lb = [];
                ub = [];
                weight = [];
                grad = [];
                hessian = [];
                tolerance = [];
                map = [];
                if any(target_I)
                    if opts.LB
                        lb = opts.LB(target_I);
                    end
                    if opts.UB
                        ub = opts.UB(target_I);
                    end
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

                [cd_comp_sorted, d_sorted] = processComp(target_comp, cd_comp, n_max, weight, grad, hessian, map, tolerance, lb, ub);
                
                CD_sorted{1,i} = cd_comp_sorted;
                D_sorted{1,i} = d_sorted;
            end
            
            if N_unique_comps == 1
                CD_sorted = CD_sorted{:};
                D_sorted = D_sorted{:};
            end
            
            function [cd_sorted,d_sorted] = processComp(target, cd, n_max, weights, grad, hessian, map, tolerance, lb, ub)
                N_comp = numel(cd);
                d = zeros(N_comp,1);
                for j = 1:N_comp
                    d(j) = distance(target, cd(j).Data, 'Mode', opts.DistanceMode, 'LB', lb, 'UB', ub, 'Weights', weights, 'Gradient', grad, 'Hessian', hessian, 'Map', map, 'Tolerance', tolerance);
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

            assert(all(obj_array(1).Component == [obj_array.Component]), 'Component Data Objects must have Common Component to make a table');
            
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

