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
        
        function t = table(obj_array)
            w=warning('off','MATLAB:structOnObject');
            t=struct2table(arrayfun(@struct, obj_array));
            warning(w);
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

