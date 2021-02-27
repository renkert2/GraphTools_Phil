classdef DynamicTypes < uint8
    % DynamicTypes is an enumeration class that defines how to cacluate a 
    % graph model's state dynamics. These class helps is used to conver a
    % graph model to a "modified graph model" as described in Section 3 of
    % "Hierarchical Model-Based Predictive Controller for a Hybrid UAV
    % Powertrain" By C.T. Aksland and A.G. Alleyne
    % Enumerated value options are:
    % (1) Energy Flow
    % (2) State Flow
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Contributors: Christopher T. Aksland and Phil Renkert
    % Association: University of Illionis at Urbana-Champaign
    % Contact: aksland2@illinois.edu and renkert2@illinois.edu
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % potential improvements:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    enumeration
        EnergyFlow (1)
        StateFlow (2)
    end
end