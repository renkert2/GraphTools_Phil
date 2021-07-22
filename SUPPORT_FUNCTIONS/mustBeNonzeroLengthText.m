function mustBeNonzeroLengthText(text)
%mustBeNonzeroLengthText Validate that value has non-zero strlength
%   mustBeNonzeroLengthText(TEXT) throws an error if TEXT does not have
%   least one character in each element.
%
%   See also mustBeText, strlength, isstring, ischar, validators.
    
% Copyright 2020 The MathWorks, Inc.

    if ~istext(text) 
        throwAsCaller(MException("MATLAB:validators:mustBeNonzeroLengthText", message("MATLAB:validators:mustBeText")));
    elseif ~all(strlength(text)>0,'all')
        throwAsCaller(MException("MATLAB:validators:mustBeNonzeroLengthText", message("MATLAB:validators:nonzeroLengthText")));
    end
end