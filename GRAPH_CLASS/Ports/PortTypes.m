classdef PortTypes < uint8
    % PortTypes is an enumeration class that defines types of component
    % connections. Components can be connected along an edge (Type_1) or
    % at a vertex (Type_2). See Figure 2.19 of Christopher T. Aksland for
    % more details (note that the Type_1 and Type_2 definitions may be the 
    % opposite in the thesis
    % Objects of this class do not need to be accessed by a user.
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Contributors: Christopher T. Aksland and Phil Renkert
    % Association: University of Illionis at Urbana-Champaign
    % Contact: aksland2@illinois.edu and renkert2@illinois.edu
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % potential improvements:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    enumeration
        Type_1 (1) % Edge Connection Type
        Type_2 (2) % Vertex Connection Type
    end
end
