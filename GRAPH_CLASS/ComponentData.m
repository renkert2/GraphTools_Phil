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
        
        function [tbl,param_table] = table(obj_array)
            pdat = obj_array(1).Data;
            param_table = table(pdat);
            pvars = string(param_table.Properties.VariableNames);
            pvi = pvars ~= "Value";
            param_table = param_table(:,pvi);
            
            params = [pdat.Sym];
            comp_fields = ["Make", "Model", "SKU", "Description"];
            varnames = [comp_fields, params];
            vartypes = [repmat("string",1,4) repmat("double",1,numel(params))];
            tbl = table('Size', [numel(obj_array), numel(varnames)], 'VariableTypes', vartypes, 'VariableNames', varnames);
            for i = 1:numel(obj_array)
                for f = comp_fields
                    val = obj_array(i).(f);
                    if isempty(val)
                        val = "";
                    end
                    tbl(i,:).(f) = val;
                end
                for j = 1:numel(params)
                    f = params(j);
                    val = obj_array(i).Data(j).Value;
                    if isempty(val)
                        val = NaN;
                    end
                    tbl(i,:).(f) = val;
                end
            end
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

