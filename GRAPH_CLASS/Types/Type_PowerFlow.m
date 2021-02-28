classdef Type_PowerFlow < Type
    % Type_PowerFlow (subclass of see Type) is used to define the 
    % expression to calculate constant coefficient powerflows in the 
    % Graph Modeling Toolbox. The powerflows must be defined in terms of 
    % the tail state 'xt', head state 'xh', and inputs 'u1','u2',...,'uN'.
    % Instatiate an object in the following form:
    % 
    % T = Type_PowerFlow(type) where type defines the expression used to
    % calulacte the power flow of an edge in a Graph object. type must
    % only be a function of the tail state xt, head state xh, and inputs 
    % 'u1','u2',...,'uN', (as a string or symbolic variable).
    % Ex: Suppose the powerflow of a vertex is calculated as c1*(u1*xt+u2*xh)
    % where c1 is a constant coefficient, xt is the tail state, xh is the
    % head state, and u1 and u2 are inputs
    % String definition: 
    % T = Type_Capacitance('u1*xt+u2*xh');
    % Symbolic definition:
    % xt = sym('xt')
    % xh = sym('xh')
    % u1 = sym('u1')
    % u2 = sym('u2')
    % T = Type_Capacitance(u1*xt+u2*xh)
      
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Contributors: Christopher T. Aksland and Phil Renkert
    % Association: University of Illionis at Urbana-Champaign
    % Contact: aksland2@illinois.edu and renkert2@illinois.edu
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % potential improvements:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
    
    properties (SetAccess = private)
        num_inputs {mustBeInteger, mustBeNonnegative} = 0 % Number of inputs incident on the power flow
    end
    
    methods
        function obj = Type_PowerFlow(type)
            [type, params_temp, num_inputs_temp] = Type_PowerFlow.parsePowerFlowVars(type); % Parses type, returns params argument structure and the number of inputs
            obj = obj@Type(type, params_temp);
            
            obj.num_inputs = num_inputs_temp;
        end
        
        
        function val = calcVal(obj, xt, xh, u) 
            % Substitutes numerical values for xt, xh, [u1...un].
            % obj can be a single object or object array of Type_PowerFlows
            
            num_args = nargin-1; % Matlab counts obj as an input arguments
            num_inputs_array = [obj.num_inputs];
            
            if num_args == 3 && ~isempty(u) % Input arguments 'u' specified
                assert(iscolumn(xt), 'xt must be column vector')
                assert(size(xt,1) == size(xh,1) &&  size(xt,1) == size(u,1), 'Dimension 1 of xt, xh, and u must be equivalent')
                if any(num_inputs_array) % if any of the Type_PowerFlows in obj require input arguments
                    assert(length(u) >= max(num_inputs_array), 'Invalid number of inputs. Number of columns of u must be greater than or equal to the maximum number of inputs');
                    val = calcVal@Type(obj, xt, xh, u);
                else % input arguments not required
                    val = calcVal@Type(obj,xt,xh);
                end
            elseif num_args == 2 || isempty(u) % input arguments 'u' not specified
                assert(iscolumn(xt), 'xt must be column vector')
                assert(size(xt,1) == size(xh,1), 'Dimension 1 of xt and xh must be equivalent')
                assert(all(num_inputs_array == 0), "Input argument 'u' required'");
                val = calcVal@Type(obj, xt, xh);
            else
                error("Invalid arguments")
            end
        end
        
        
        function type = updateInputs(obj, type)
            [type, params_temp, num_inputs_temp] = Type_PowerFlow.parsePowerFlowVars(type);
            obj.params = params_temp;
            obj.num_inputs = num_inputs_temp;
        end
        
        function init(obj, type)
            type = obj.updateInputs(type);
            init@Type(obj, type);
        end
        
    
end

methods(Static)
    function [type, params, num_inputs] = parsePowerFlowVars(type)
        % Returns symbolic variables, params argument structure, and number of inputs incident on a type
        
        [type,num_inputs] = countInputVars(type);
        params = calcVars(num_inputs);
        
        function [type, num_inputs] = countInputVars(type)
            % Converts type to a string and uses regular expression to count the number of inputs required
            
            if isa(type, 'sym') || isa(type, 'char')
                type = string(type);
            elseif ~isa(type,'string')
                error("PowerFlow Type must be of type string, char, or sym")
            end
            
            patt_single_input = 'u(?!\d)'; % Finds 'u' occurences in expression
            type = regexprep(type, patt_single_input, "u1"); % Replace occurences of 'u' with 'u1'

            patt_list = '(?<=u)\d+';
            input_list = double(regexp(type, patt_list, 'match')); % Finds 'u[12...N]' occurences in expression
            input_list_flag = ~isempty(input_list); % Flag indicating if any inputs were found
            
            if input_list_flag
                num_inputs = max(input_list);
            else
                num_inputs = 0;
            end
        end
        
        function [params] = calcVars(num_inputs)
            state_vars = [sym('xt'), sym('xh')];
            if num_inputs >= 1
                input_vars =sym('u', [1, max(num_inputs,2)]);
                params = {state_vars(1), state_vars(2), input_vars};
            elseif num_inputs == 0
                params = {state_vars(1), state_vars(2)};
            end
        end
    end
    
end
end

