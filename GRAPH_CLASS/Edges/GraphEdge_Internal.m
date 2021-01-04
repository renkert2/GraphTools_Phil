classdef GraphEdge_Internal < GraphEdge
    %GraphEdge_Internal Summary of this class goes here
    %   Detailed explanation goes here
    properties
        PowerFlow Type_PowerFlow = Type_PowerFlow.empty()
        Input GraphInput = GraphInput.empty()
        Port
        Coefficient (:,1) double = 0
        TailVertex (1,1) GraphVertex = GraphVertex()
    end
    
    methods
        function obj = GraphEdge_Internal(varargin) % constructor method
            obj@GraphEdge(varargin{:}); % calls the superclass constructor
        end
        
        function set.TailVertex(obj, tail_vertex)
            assert(numel(tail_vertex)<=1, 'Edge can have only one tail vertex');
            obj.TailVertex = tail_vertex;
        end
        
                
        function set.PowerFlow(obj, power_flow) % Ensure PowerFlows are column vector
            if ~iscolumn(power_flow)
                power_flow = power_flow';
            end
            obj.PowerFlow = power_flow;
        end
        
        function set.Input(obj, input) % Ensure inputs are column vector
            if ~iscolumn(input)
                input = input';
            end
            obj.Input = input;
        end
                 
        function set.Port(obj, port) % Ensure ports are column vector
            if ~iscolumn(port)
                port = port';
            end
            obj.Port = port;
        end
        
                
        function set.Coefficient(obj, coeff) % Ensure coefficients are column vector
            if ~iscolumn(coeff)
                coeff = coeff';
            end
            obj.Coefficient = coeff;
        end
    end
end
