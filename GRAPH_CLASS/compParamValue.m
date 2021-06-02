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
        
        function [obj_filtered, I] = filterComponent(obj_array, component)
            I = vertcat(obj_array.Component) == component;
           obj_filtered = obj_array(I); 
        end
        
        function obj_filtered = filterSym(obj_array, sym)
            obj_filtered = obj_array([obj_array.Sym] == sym); 
        end
        
        function [obj_pairs, I] = getPairs(obj_array_1, obj_array_2)
            % Each row of obj_pairs is a parameter pair
            % obj_pairs(1,:) = obj_array_1(I(1,:))
            % obj_pairs(2,:) = obj_array_2(I(2,:))
            obj_pairs = compParamValue.empty(0,2);
            I = double.empty(0,2);
            pair_cnt = 0;
            for i = 1:numel(obj_array_1)
                obj_1 = obj_array_1(i);
                f_sym = obj_1.Sym == [obj_array_2.Sym];
                f_comp = obj_1.Component == [obj_array_2.Component];
                f_combined = f_sym & f_comp;
                n = sum(f_combined);
                if n == 1
                    pair_cnt = pair_cnt + 1;
                    I(pair_cnt, :) = [i, find(f_combined)];
                    obj_pairs(pair_cnt, :) = [obj_1, obj_array_2(f_combined)];
                end
            end
        end
        
        function d = distance(target, candidate, opts)
            % Optionally pass weights as a struct with key = sym, value =
            % weight
            arguments 
                target compParamValue
                candidate compParamValue
                opts.Mode string = "Norm" % Modes: Norm, WeightedNorm, TaylorSeries
                opts.Weights = []
                opts.Gradient double = [] % Must be in same order as target
                opts.Hessian double = [] % Must be in same order as target
            end
            
            [pairs, I_pairs] = getPairs(target, candidate); % Nx2 compParam array
            
            target_vals = vertcat(pairs(:,1).Value);
            candidate_vals = vertcat(pairs(:,2).Value);
            
            switch opts.Mode
                case "Norm"
                    normalized_distance = (candidate_vals./target_vals - 1);
                    d = norm(normalized_distance);
                case "WeightedNorm"
                    if isa(opts.Weights, 'struct')
                        weights_struct = opts.Weights;
                        weights = ones(size(target_vals));
                        syms = [pairs(:,1).Sym];
                        for i = 1:numel(syms)
                            if isfield(weights_struct, syms(i))
                                weights(i) = weights_struct.(syms(i));
                            end
                        end
                    else
                        weights = opts.Weights(I_pairs(:,1));
                    end
                    
                    normalized_distance = abs(weights).*(candidate_vals./target_vals - 1);
                    d = norm(normalized_distance);
                case "TaylorSeries"
                    dx = candidate_vals - target_vals;
                    I = I_pairs(:,1);
                    
                    if ~isempty(opts.Gradient)
                        G = opts.Gradient(I);
                    else
                        G = zeros(numel(I),1);
                    end
                    
                    if ~isempty(opts.Hessian)
                        H = opts.Hessian(I,I);
                    else
                        H = zeros(numel(I));
                    end
                    
                    df = G'*dx + (1/2)*dx'*H*dx; % Approximate change in function
                    d = df;
                otherwise
                    error("Invalid Mode.  Must be Norm, WeightedNorm, or TaylorSeries");
            end

        end
        
        function [obj_array_sorted, I] = sort(obj_array, varargin)
            [~,I] = sort([obj_array.Sym]);
            obj_array_sorted = obj_array(I);
        end
        
        function s = makeValStruct(obj_array, X)
            if nargin == 1
                X = vertcat(obj_array.Value);
            end
            
            s = struct();
            for i = 1:numel(obj_array)
                obj = obj_array(i);
                s.(obj.Component).(obj.Sym) = X(i);
            end
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

