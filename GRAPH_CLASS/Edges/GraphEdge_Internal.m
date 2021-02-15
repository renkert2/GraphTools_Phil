classdef GraphEdge_Internal < GraphEdge
    %GraphEdge_Internal Summary of this class goes here
    %   Detailed explanation goes here
    properties
        PowerFlow Type_PowerFlow = Type_PowerFlow.empty()
        Input GraphInput = GraphInput.empty()
        Coefficient (:,1) = 0
        TailVertex GraphVertex = GraphVertex.empty()
        HeadVertex GraphVertex = GraphVertex.empty()
    end
    
    methods
        function obj = GraphEdge_Internal(varargin) % constructor method
            obj@GraphEdge(varargin{:}); % calls the superclass constructor
        end
        
        function set.TailVertex(obj, tail_vertex)
            assert(numel(tail_vertex)<=1, 'Edge can have only one tail vertex');
%             if ~isempty(obj.HeadVertex)
%                 if ~isConnectable(obj.HeadVertex.VertexType, tail_vertex.VertexType)
%                     warning('Tail vertex is not connectable with existing head vertex')
%                 end
%             end
            obj.TailVertex = tail_vertex;
        end
        
        function set.HeadVertex(obj, head_vertex) % Overrides set.HeadVertex of GraphEdge
            assert(numel(head_vertex)<=1, 'Edge can have only one head vertex');
%             if ~isempty(obj.TailVertex)
%                 if ~isConnectable(obj.TailVertex.VertexType, head_vertex.VertexType)
%                     warning('Head vertex is not connectable with existing tail vertex')
%                 end
%             end
            obj.HeadVertex = head_vertex;
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
                          
        function set.Coefficient(obj, coeff) % Ensure coefficients are column vector
            if ~iscolumn(coeff)
                coeff = coeff';
            end
            obj.Coefficient = coeff;
        end
    end
end
