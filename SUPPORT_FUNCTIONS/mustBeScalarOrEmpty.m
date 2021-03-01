function mustBeScalarOrEmpty(A)
%MUSTBESCALAROREMPTY Validate that value is a scalar or is empty
%   MUSTBESCALAROREMPTY(A) throws an error if A is not a scalar or is not
%   empty. MATLAB calls isscalar to determine if A is a scalar, it calls
%   isempty to determine if A is empty.
%
%   Class support:
%   All MATLAB classes
%
%   See also: ISSCALAR, ISEMPTY.
        
%   Copyright 2020 The MathWorks, Inc.

    if ~isscalar(A) && ~isempty(A)
        throwAsCaller(MException(message('MATLAB:validators:mustBeScalarOrEmpty')));
    end
end
