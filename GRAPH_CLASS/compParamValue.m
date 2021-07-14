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
                opts.LB double = []
                opts.UB double = []
                opts.Weights = []
                opts.Gradient double = [] % Must be in same order as target
                opts.Hessian double = [] % Must be in same order as target
                opts.Map double = [] % Maps target parameter values to their respective location in the Gradient/Hessian
                opts.Tolerance = Inf % Determines how far away the second-order estimation can look.  Each normalized error must be less than the tolerance.
            end
            
            N_target = numel(target);
            sz_target = size(target);
            
            [pairs, I_pairs] = getPairs(target, candidate); % Nx2 compParam array
            
            
            target_vals = vertcat(pairs(:,1).Value);
            candidate_vals = vertcat(pairs(:,2).Value);
            
            sz_pairs = size(target_vals);
            
            % Parse Upper and Lower Bound Options
            if opts.LB
                if isscalar(opts.LB)
                    LB = repmat(opts.LB, sz_target);
                else
                    assert(numel(opts.LB) == N_target, "Number of vars in LB must match size of target");
                    LB = opts.LB;
                end
            else
                LB = -Inf(sz_target);
            end       
            lb = LB(I_pairs(:,1));
            
            if opts.UB
                if isscalar(opts.UB)
                    UB = repmat(opts.UB, sz_target);
                else
                    assert(numel(opts.UB) == N_target, "Number of vars in UB must match size of target");
                    UB = opts.UB;
                end
            else
                UB = Inf(sz_target);
            end
            ub = UB(I_pairs(:,1));
            
            in_bounds = (candidate_vals >= lb) & (candidate_vals <= ub);
            if ~all(in_bounds)
                d = NaN;
                return
            end
            
            switch opts.Mode
                case "Norm"
                    normalized_distance = (candidate_vals./target_vals - 1);
                    d = norm(normalized_distance);
                case "WeightedNorm"
                    if isa(opts.Weights, 'struct')
                        weights_struct = opts.Weights;
                        weights = ones(sz_pairs);
                        syms = [pairs(:,1).Sym];
                        for i = 1:numel(syms)
                            if isfield(weights_struct, syms(i))
                                weights(i) = weights_struct.(syms(i));
                            end
                        end
                    else
                        assert(numel(opts.Weights) == N_target, "Size of weights must match number of targets");
                        weights = opts.Weights;
                        weights = weights(I_pairs(:,1));
                    end

                    normalized_distance = abs(weights).*(candidate_vals./target_vals - 1);
                    d = norm(normalized_distance);
                case "TaylorSeries"
                    if isa(opts.Tolerance, 'struct')
                        tol_struct = opts.Tolerance;
                        tols = inf(sz_pairs);
                        syms = [pairs(:,1).Sym];
                        for i = 1:numel(syms)
                            if isfield(tol_struct, syms(i))
                                tols(i) = tol_struct.(syms(i));
                            end
                        end
                    else
                        if isnumeric(opts.Tolerance) && isscalar(opts.Tolerance)
                            tols = repmat(opts.Tolerance, sz_target);
                        elseif isempty(opts.Tolerance)
                            tols = inf(sz_target);
                        else
                            assert(numel(opts.Tolerance) == N_target, "Size of weights must match number of targets");
                            tols = opts.Tolerance;
                        end
                        tols = tols(I_pairs(:,1));
                    end
                    
                    
                    dx = candidate_vals - target_vals;

                    normalized_distance = (candidate_vals./target_vals - 1);
                    valid = all(abs(normalized_distance) < tols);
                    
                    if valid
                        I_x = I_pairs(:,1);

                        if ~isempty(opts.Gradient)
                            G = opts.Gradient;
                            N_G = numel(G);
                        else
                            N_G = 0;
                        end

                        if ~isempty(opts.Hessian)
                            H = opts.Hessian;
                            [N_H_1, N_H_2] = size(H);
                            assert(N_H_1 == N_H_2, "Hessian must be square");
                            N_H = N_H_1;
                        else
                            N_H = 0;
                        end

                        if all([N_G, N_H])
                            assert(N_G == N_H, 'Gradient and Hessian must be of same size');
                            N_X = max([N_G, N_H]); % Total number of inputs required to evaluate function approximation
                        elseif N_H == 0
                            N_X = N_G;
                            H = zeros(N_X);
                        elseif N_G == 0
                            N_X = N_H;
                            G = zeros(N_X,1);
                        end


                        if N_X < N_target
                            error("Size of gradient/hessian must be at least as large as the number of target variables")
                        else
                            if ~isempty(opts.Map)
                                assert(numel(opts.Map) == N_target, "Size of map must match number of targets")
                                I_X = opts.Map(I_x);
                            else
                                assert(N_X == N_target, "Map option required if gradient/hessian is larger than the number of target variables")
                                I_X = I_x;
                            end
                        end

                        dX = zeros(N_X,1);
                        dX(I_X) = dx;

                        df = G'*dX + (1/2)*dX'*H*dX; % Approximate change in function
                        d = df;
                    else
                        d = NaN;
                    end
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

