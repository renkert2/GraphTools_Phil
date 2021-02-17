function mustBeNonnegativeOrSym(p)
    if ~isa(p,"sym")
        mustBeNonnegative(p);
    end
end