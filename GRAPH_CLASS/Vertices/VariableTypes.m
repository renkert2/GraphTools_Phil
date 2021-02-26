classdef VariableTypes < uint8
    % VariableTypes is an enumeration class that distinguishes between the 
    % bond graph nomenclature of effort and flow variables in the Graph 
    % Modeling Toolbox. Although graph models and bond graphs differ, the
    % concept of effor and flow variables are useful in understand graph
    % construction compatability.
    % Enumerated value options are:
    % (0) Abstract
    % (1) Effort
    % (2) Flow
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Contributors: Christopher T. Aksland and Phil Renkert
    % Association: University of Illionis at Urbana-Champaign
    % Contact: aksland2@illinois.edu and renkert2@illinois.edu
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % potential improvements:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    enumeration
        Abstract (0)
        Effort (1)
        Flow (2)
    end
end