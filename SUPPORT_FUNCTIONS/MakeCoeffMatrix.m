function [Coeff,TypeSys] = MakeCoeffMatrix (All,TypeAll,numType)
    % MakeCoeffMatrix is a function develops coefficient matrices for the
    % component graphs. Information stored in as Type and Value is
    % converted to matrices that are used the system generation code
    
    %%% INPUTS
    % Comp - Cell Array of Comp_Graph classes
    % ic - vector map of unique types to all types. This is the output of 
    %      the UniqueString function
    % indicator - String prefix for the type and coefficient matrix 'C' and
    %             'P' are currently supported
    % type - cell array of types
    % val - cell array of values

    %%% OUTPUTS
    % Comp - Cell array of Comp_Graph classes with updated coefficient
    % matrices
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Author: Christopher T. Aksland
    % Association: University of Illionis at Urbana-Champaign
    % Contact: aksland2@illinois.edu
    % Revision History:
    % 7/29/2020 - Function creation
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    [~, IB, IC] = unique(vertcat(TypeAll(:).Val_Sym)); % get unique capactitance types
    Coeff = zeros(length(All),max(IC)); % initialize stacked capacitance coefficient matrix
    if any(arrayfun(@(x) isa(x.Coefficient, 'sym'), All))
        Coeff = sym(Coeff);
    end
    
    j = 0;
    for i = 1:length(All) % loop through all vertices and assign coefficients
        Coeff(i,IC(i+j:i-1+numType(i)+j)) = All(i).Coefficient; % assign values
        j = j+numType(i)-1; % j accounts for instances where a vertex may have multiple coefficients
    end
    
    
    
    TypeSys = copy(TypeAll(IB)); % get the system capacitance types, use "copy" to create unique instance/handle

    
end