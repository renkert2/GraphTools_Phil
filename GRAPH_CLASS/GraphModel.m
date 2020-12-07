classdef GraphModel < matlab.mixin.Copyable
    %GRAPHMODEL Contains all of the operations useful for working with the
    %graph models.  Most of the code will go here.  
    %   Detailed explanation goes here
    properties
        Graph Graph = Graph.empty()
    end
    
    methods
        function obj = GraphModel(varargin)
            if nargin == 1
                obj.Graph = varargin{1};
            end
        end
        
        function plot(obj)
            % basic digraph plotting.
%             figure
            G = digraph(obj.Graph.E(:,1),obj.Graph.E(:,2));
            h = plot(G);
            labeledge(h,obj.Graph.E(:,1)',obj.Graph.E(:,2)',1:size(obj.Graph.E,1));
        end
        
        function Modify(obj)
        
        end
        
        function Simulate(obj)
            
        end
        
        function init(obj)
            % placeholder
        end
    end
    
    
    
    methods(Static)
        function g_sys = Combine(G, ConnectV, ConnectE) % Create a new GraphModel Object
            % Algorithm 1
            % Algorithm 2
        end
        
    end
        
end

