classdef CalcPMethods
    enumeration
        Default % vectorized calculation of every power flow
        Edges % calculates powerflows independently for each edge.  Avoids some numerical issues with complicated powerflows. 
    end
end