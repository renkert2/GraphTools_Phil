classdef Graph < matlab.mixin.Copyable
    %GRAPH Fundamental Graph Class
    %   Keep pure and simple: A graph consists of Vertices, Edges, and the
    %   Connection Matrix / Incidence Matrix
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Author: Christopher T. Aksland and Phil Renkert
    % Association: University of Illionis at Urbana-Champaign
    % Contact: aksland2@illinois.edu
    % Revision History:
    % 12/5/2020 - Consider added Matlab "digraph" or adjcancy matrix as a 
    %           propertey.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % potential improvements:
    % Consider added Matlab "digraph" or adjcancy matrix as a propertey.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    properties
        Vertices GraphVertex_Super
        Edges GraphEdge_Super
        E % Graph Edge Matrix
    end
    
    properties (SetAccess = private)
        % Incidence Matrix
        M 
        
        % Number of vertices
        Nv
        
        % Number of edges
        Ne
        
        % Number of inputs
        Nu
        
        % Number of extternal vertices
        Nev
        
        % Number of external edges
        Nee
      
    end
    
    methods
        function plot(obj)
            % basic digraph plotting.
            figure
            G = digraph(obj.E(:,1),obj.E(:,2));
            h = plot(G);
            labeledge(h,obj.E(:,1)',obj.E(:,2)',1:size(obj.E,1));               
        end
    end
end

