classdef compParamValue
    % Value class used by compParam to load parameter values representing
    % physical components.  
    
    properties
        Sym string
        Value double
        Unit string
        Description string
        
        Component string % Name of component object corresponding to this parameter
        % May want to add ComponentData parent here later
    end
    
    methods
        function obj = compParamValue(sym, varargin)
            if nargin == 1
                obj.Sym = sym;
            elseif nargin > 1
                obj.Sym = sym;
                obj = my_inputparser(obj, varargin{:});
            end
        end
        
        function obj_filtered = filterComponent(obj_array, component)
           obj_filtered = obj_array([obj_array.Component] == component); 
        end
        
        function obj_pairs = getPairs(obj_array_1, obj_array_2)
            % Each row of obj_pairs is a parameter pair
            obj_pairs = compParamValue.empty(0,2);
            pair_cnt = 0;
            for i = 1:numel(obj_array_1)
                obj_1 = obj_array_1(i);
                f_sym = obj_1.Sym == [obj_array_2.Sym];
                f_comp = obj_1.Component == [obj_array_2.Component];
                f_combined = f_sym & f_comp;
                n = sum(f_combined);
                if n == 1
                    pair_cnt = pair_cnt + 1;
                    obj_pairs(pair_cnt, :) = [obj_1, obj_array_2(f_combined)];
                end
            end
        end
        
        function d = distance(target, candidate)
            pairs = getPairs(target, candidate); % Nx2 compParam array
            
            target_vals = vertcat(pairs(:,1).Value);
            candidate_vals = vertcat(pairs(:,2).Value);
            
            normalized_distance = candidate_vals./target_vals - 1;
            d = norm(normalized_distance);
        end
        
        function t = table(obj_array)
            w=warning('off','MATLAB:structOnObject');
            t=struct2table(arrayfun(@struct, obj_array));
            warning(w);
        end
        
        function dispTable(obj_array)
            % Modifies method from Mixin.Custom Display
            tbl = table(vertcat(obj_array.Sym), vertcat(obj_array.Value), vertcat(obj_array.Unit), vertcat(obj_array.Description), vertcat(obj_array.Component),...
                'VariableNames', ["Sym", "Value", "Unit", "Description", "Component"]);
            disp(tbl);
        end

    end

    methods (Static)
        function obj_array = importFromTable(param_table)
            % Construct array of compParamValues from table, whos
            % VariableNames correspond to compParamValue properties.
            % First column must be Sym
            
            N = height(param_table);
            propnames = string(param_table.Properties.VariableNames(2:end));
            syms = param_table{:,1};
            
            obj_array = compParamValue.empty(N,0);
            for i = 1:N
                obj = compParamValue(syms(i));
                for j = 1:numel(propnames)
                    field = propnames(j);
                    obj.(field) = param_table(i,:).(field);
                end
                obj_array(i,1) = obj;
            end
        end
        
        function obj_array = importFromStruct(param_struct)
            % Construct array of compParamValues from struct array or cell array 
            % of structs
            
            N = numel(param_struct);  
            obj_array = compParamValue.empty(N,0);
            for i = 1:N
                obj = compParamValue();
                if isa(param_struct, 'cell')
                    ps = param_struct{i};
                else
                    ps = param_struct(i);
                end
                field_names = string(fields(ps))';
                for f = field_names
                    obj.(f) = ps.(f);
                end
                obj_array(i,1) = obj;
            end
        end
    end
end

