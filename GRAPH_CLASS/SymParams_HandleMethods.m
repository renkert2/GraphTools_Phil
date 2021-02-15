classdef SymParams_HandleMethods
    enumeration
        AugmentMatlabFunctions % Handles symbolic parameters by adding an extra input to CalcX functions
        SubstituteNumericalValues % Handles symbolic parameters by compiling CalcX functions with values in Model.SymParams_Vals
    end
end