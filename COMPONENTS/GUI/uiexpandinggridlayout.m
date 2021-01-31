function gridComponent = uiexpandinggridlayout(varargin)
%UIGRIDLAYOUT Create a grid layout container
%   g = UIGRIDLAYOUT creates a grid layout in a new figure and returns the
%   GridLayout object. MATLAB calls the uifigure function to create the
%   figure.
%
%   g = UIGRIDLAYOUT(parent) creates the grid layout in the specified
%   parent container. The parent container can be a figure created using
%   the uifigure function, or one of its child containers such as a panel.
%
%   g = UIGRIDLAYOUT(sz) creates the grid layout, and specifies the size as
%   a two-element vector sz. The first element of sz is the number of rows
%   in the grid. The second element is the number of columns.
%
%   g = UIGRIDLAYOUT(parent, sz) creates the grid layout in the specified
%   parent container with the specified size.
%
%   g = UIGRIDLAYOUT( ___ ,Name,Value) specifies grid layout properties
%   using one or more Name,Value pair arguments. Use this option with any
%   of the input argument combinations in the previous syntaxes.
%
%   Example 1: Create a Grid Layout
%      % Create a grid layout
%      gridlayout = uigridlayout;
%
%   Example 2: Specify Parent Object for Grid Layout
%      % Specify a UI figure as the parent object for a grid layout.
%      fig = uifigure;
%      gridlayout = uigridlayout(fig);
%
%   Example 3: Specify Size for Grid Layout
%      % Create a grid layout with 4 rows and 5 columns.
%      gridlayout = uigridlayout([4 5]);
%
%   See also UIFIGURE, UIPANEL, UITAB

%   Copyright 2018 The MathWorks, Inc.

% Check first two input arguments for grid size

className = 'ExpandingGrid';

messageCatalogID = 'uiexpandinggridlayout';

try
    updatedInput = handleSizeInput(varargin);
catch ex
    error('MATLAB:ui:GridLayout:unknownInput', ...
        ex.message);
end

try
    gridComponent = matlab.ui.control.internal.model.ComponentCreation.createComponent(...
        className, ...
        messageCatalogID,...
        updatedInput{:});
catch ex
    error('MATLAB:ui:GridLayout:unknownInput', ...
        ex.message);
end

    function input = handleSizeInput(input)
        % handleSizeInput returns an updated input if the first or second input is
        % a valid grid size.
        
        
        if isempty(input) || isPossiblePropertyName(input{1})
            % Input has neither parent nor size
            return;
        end
        
        if isPossibleParent(input{1})
            % First argument might be valid parent, check second input
            % (Note: more detailed checks for parent done during object
            % creation)
            
            if length(input) > 1 && ~isPossiblePropertyName(input{2})
                % Second input is a size
                if isValidGridSize(input{2})
                    input = updateInput(input, 2);
                    return;
                else
                    error(message('MATLAB:ui:containers:mustBe1x2RowOfPositiveNumbers'));
                end
                
            else
                % Second input is not a size
                return;
            end
            
        end
        
        if isPossibleSize(input{1})
            % First argument was intended as a size
            
            if isValidGridSize(input{1})
                input = updateInput(input, 1);
                return;
            else
                error(message('MATLAB:ui:containers:mustBe1x2RowOfPositiveNumbers'));
            end
        end
        
        if isobject(input{1}) && ~isvalid(input{1})
            % First argument is a deleted handle
            error(message('MATLAB:ui:components:invalidObject', 'Parent'));
        end
        
    end

    function tf = isPossibleParent(input)
        % isPossibleParent returns true if input might be a valid parent. More
        % detailed checks are done during object creation.
        tf = (isobject(input) && isgraphics(input));
    end

    function tf = isPossiblePropertyName(x)
        % isPossibleParameter Name returns true if input is a char or string
        tf = isstring(x) || ischar(x);
    end

    function tf = isPossibleSize(x)
        % isPossibleSize returns whether the input is numeric or a cell
        % with numeric values
        
        isNumericCell = iscell(x) && all(cellfun(@(k)isnumeric(k), x));
        tf = isnumeric(x) || isNumericCell;
    end

    function tf = isValidGridSize(x)
        % validGridSize returns whether input is a valid grid size.
        
        try
            validateattributes(x, {'numeric'}, ...
                {'row', 'numel', 2, 'real', 'finite', 'nonnan', ...
                'positive', 'integer'});
            tf = true;
        catch
            tf = false;
        end
    end

    function input = updateInput(input, index)
        % updateInput converts grid size input to Name/Value pairs
        x = input{index};
        
        nRows = repmat({'1x'}, 1, x(1));
        nColumns = repmat({'1x'}, 1, x(2));
        
        input = [input(1:index-1) ...
            {'RowHeight', nRows, 'ColumnWidth', nColumns} ...
            input(index+1:end)];
    end
end



