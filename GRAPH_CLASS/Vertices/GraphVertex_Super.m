classdef GraphVertex_Super < matlab.mixin.Heterogeneous
    %GRAPHVERTEX_Super Super Class for all vertex types (i.e. State,
    %External)
    %   The matlab.mixin.Heterogeneous type will allow us to create heterogeneous object
    %   arrays of all the child classes.  
    properties
        Description
    end
end

