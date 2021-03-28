function mustBeNumericOrSym(p)
    if ~isa(p,"sym")
        mustBeNumeric(p);
    end
end