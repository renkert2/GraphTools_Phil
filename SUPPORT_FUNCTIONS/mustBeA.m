function mustBeA(A,C)
%MUSTBEA Validate that value comes from one of the specified classes
%   MUSTBEA(A,C) compares A with a list of class names in C and throws an
%   error if the class of A is not one of the classes or a subclass of one
%   of the classes. C can be a string array, a character vector, or a
%   cell array of character vectors.
%
%   MATLAB calls the isa to determine if A is an object of a class.
%
%   Class support:
%   All MATLAB classes
%
%   See also: ISA.
        
%   Copyright 2020 The MathWorks, Inc.

    arguments
        A
        C {mustBeNonzeroLengthText}
    end
    
    % empty C, including {} and empty string, makes any value valid
    if isempty(C)
        return;
    end
    
    C = string(C);
    
    if ~any(arrayfun(@(cls)isa(A,cls), C), 'all')
        if numel(C) <= 6
            msgIDs = [...
                "MATLAB:validators:OneType",...
                "MATLAB:validators:TwoTypes",...
                "MATLAB:validators:ThreeTypes",...
                "MATLAB:validators:FourTypes",...
                "MATLAB:validators:FiveTypes",...
                "MATLAB:validators:SixTypes"...
                ];
            messageObject = message(msgIDs(numel(C)), C{1:end});
            throwAsCaller(MException(message('MATLAB:validators:mustBeA', messageObject.getString)));
        else
            formatedList = createPrintableList(C);
            throwAsCaller(MException(message('MATLAB:validators:mustBeA', formatedList)));
        end
    end
end


% LocalWords:  validators
