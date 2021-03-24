classdef SymParams < matlab.mixin.Copyable
    %SYMPARAMS Contains list of SymParams Symbolic Variables, Default Values, and the number of SymParams
    % Note that this is an immutable Value class
    properties (SetAccess = private)
        Syms (:,1) sym = sym.empty()
        Vals (:,1) double = []
        N (1,1) double = 0
    end
    
    methods
        function obj = SymParams(varargin)
            if nargin == 0
                % Do Nothing
            elseif nargin == 1
                arg = varargin{1};
                if isa(arg, 'cell')
                    syms = sym.empty();
                    vals = [];
                    
                    for i = 1:numel(arg) % Loop over all cell elements
                        elem = arg{i};
                        if isa(elem, 'symParam') % For each symParam property,
                            syms(end+1,1) = elem; % Add the symParam to sym_params list
                            vals(end+1,1) = double(elem); % Add the default value of symParam to sym_params_vals list
                        end
                    end
                    
                    obj.Syms = syms;
                    obj.Vals = vals;
                    obj.N = numel(vals);
                    
                elseif isa(arg, 'sym')
                    obj.Syms = arg;
                    obj.N = numel(arg);
                    obj.Vals = zeros(obj.N, 1);
                    
                else
                    error("Single argument call to SymParams constructor must be cell or sym array");
                end
                
            elseif nargin == 2
                n = numel(varargin{1});
                assert(numel(varargin{2}) == n, "Vals must be of length %d", n);
                obj.Syms = varargin{1};
                obj.Vals = varargin{2};
                obj.N = n;
   
            else
                error("SymParams constructor requires single cell array of symParams, single sym array (Vals default to zero), or two arguments specifying syms and vals")
            end   
        end
        
        function sym_params = join(obj_array)
            syms = vertcat(obj_array.Syms);
            [syms, ia, ~] = unique(syms);
            
            vals = vertcat(obj_array.Vals);
            vals = vals(ia);
            
            sym_params = SymParams(syms,vals);
        end
        
        function append(obj, sym_param)
            assert(isa(sym_param, 'symParam'), 'Argument must be a symParam object');

            obj.Vals = [obj.Vals; double(sym_param)]; % Add the default value of symParam to sym_params_vals list
            obj.Syms = [obj.Syms; sym_param]; % Add the symParam to sym_params list
            obj.N = obj.N+1;
        end
        
        function prepend(obj, sym_param)
            assert(isa(sym_param, 'symParam'), 'Argument must be a symParam object');            
            
            obj.Vals = [double(sym_param); obj.Vals];
            obj.Syms = [sym_param; obj.Syms];
            obj.N = obj.N+1;
        end
        
        function l = isempty(obj)
            l =  builtin('isempty', obj);
            if ~l 
                l = obj.N == 0;
            end
        end
        
        function setVal(obj, param, val)
            % setVal sets the value corresponding to symParam 'param' to 'val'
            % - param: sym, char, or string identifying the symbolic variable
            % - val: double, replaces element in obj.Vals corresponding to the param
            
            if isa(param, 'char') || isa(param, 'string')
                param = sym(param);
            end
            
            i = arrayfun(@(p) isequal(p, param), obj.Syms);
            if any(i)
                obj.Vals(i) = val;
            else
                error("%s not in SymParams Syms vector", string(param));
            end
        end
    end
end

