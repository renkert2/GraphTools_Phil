classdef Graph < matlab.mixin.Copyable
    %GRAPH Fundamental Graph Class
    %   Keep pure and simple: A graph consists of Vertices, Edges, and the
    %   Connection Matrix / Incidence Matrix
    
    properties
        Vertices Graph_Vertex
        Edges Graph_Edge
        E % Graph Edge Matrix
    end
    
    properties (SetAccess = private)
        M % Incidence Matrix
        
        Ne
        Nv
        ...
    end
    
    methods
        function plot(obj)
        end
    end
end

