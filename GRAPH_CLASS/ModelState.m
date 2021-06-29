classdef ModelState
    %MODELSTATE Value class which contains system state values at an
    %instance in time.  Useful for subclassing.
    properties
        x double % Dynamic state values
        u double % Input values
        d double % Disturbance values
        y double % Output values
    end
end

