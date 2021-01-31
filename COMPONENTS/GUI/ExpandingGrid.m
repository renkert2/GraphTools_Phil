classdef ExpandingGrid < matlab.ui.container.GridLayout
    
    
    properties
%         AddButton (1,1) matlab.ui.control.Button
%         DelButton (:,1) matlab.ui.control.Button
%         Fields    (:,1) MultiField
%         Labels    (:,1) matlab.ui.control.Label
    end
    
    
    methods 
        function obj = ExpandingGrid(varargin) % constructor method
            obj@matlab.ui.container.GridLayout(varargin{:}); % calls the superclass constructor
%             obj.init();
        end
        
        function init(obj) 
%             obj.Labels(1) = uilabel(obj);
%             obj.Labels(1).HorizontalAlignment = 'center';
%             obj.Labels(1).Layout.Row = 1;
%             obj.Labels(1).Layout.Column = 2;
%             obj.Labels(1).Text = 'Property Name';
            
            
% % % %             obj.AddButton = uibutton(obj, 'push');
% % % % %             obj.AddButton.ButtonPushedFcn = createCallbackFcn(obj, @AddButtonButtonPushed, true);
% % % %             obj.AddButton.Layout.Row = 2;
% % % %             obj.AddButton.Layout.Column = 1;
% % % %             obj.AddButton.Text = 'Add Property';
            
        end
        
    end
    
    methods(Access = 'private')
        
        function implicitSize = computeImplicitSizeForSingleDirectionFromLength(obj, userSetSize, contentDrivenLength)
            
            % Number of rows or columns last set explicitly by user
            userSetLength = length(userSetSize);
            
            % Determine RowHeight or ColumnWidth
            if contentDrivenLength <= userSetLength
                % Use the user's explicitly set value
                implicitSize = userSetSize;
                
            elseif contentDrivenLength > userSetLength
                % Add implicit rows as needed to contain all components
                sizeDiff = contentDrivenLength - userSetLength;
                newSizes = repmat({'1x'},1,sizeDiff);
                implicitSize = [userSetSize, newSizes];
            end
            
            
        end
        
        function contentDrivenLength = computeContentDrivenLengthForSingleDirection(obj, direction, varargin)
            % Computes the number of implicit rows or columns. 
            % See computeImplicitGridSize
            %
            % Input: 
            %        -direction -> either 'Row' or 'Column'
            %        -varagin -> child component that is either being added or removed
            
            
            % Number of rows or columns containing all components
            % Initialize to 0 for a grid with no children
            contentDrivenLength = 0;
            
            gridChildren = obj.Children;
            if ~isempty(gridChildren)
                % Determine the smallest size that contains all the
                % children in this direction.
                % Children can span multiple cells, so use the last spanned
                % cell.
                %
                % Note: [g.Children.Layout] does not work because when
                % there is a mix of components because Layout is not
                % defined on the common base class of GBT and std/hmi
                % components, but instead mixed in each individual class.
                
                numChildren = length(gridChildren);
                location = zeros(numChildren, 1);
                
                if nargin == 3
                    childRemoved = varargin{1};
                    for k = 1:numChildren
                        % The Layout property value of the child being
                        % removed should not be accounted for 
                        % in determining the content driven size
                        if childRemoved ~= gridChildren(k)
                            layoutObj = gridChildren(k).Layout;
                            location(k) = layoutObj.(direction)(end);
                        end
                    end
                else
                
                    for k = 1:numChildren
                        layoutObj = gridChildren(k).Layout;
                        location(k) = layoutObj.(direction)(end);
                    end
                end
                contentDrivenLength = max(location);                
            end
        end
        
        % Update Content Driven Size
        function updateContentDrivenLengthForSingleDirection(obj, direction, action, varargin)
            child = [];
            childLayout = [];
            if nargin == 4
                child = varargin{1};
                childLayout = child.Layout;
            end
            
            switch action
                case 'LayoutChanged'
                    childValue = childLayout.(direction)(end);
                    if childValue > obj.ContentDrivenSize.(direction)
                        obj.ContentDrivenSize.(direction) = childValue;
                    elseif childValue < obj.ContentDrivenSize.(direction)
                        obj.ContentDrivenSize.(direction) = obj.computeContentDrivenLengthForSingleDirection(direction);
                    end
                    
                case 'childAdded'
                    childValue = childLayout.(direction)(end);
                    if childValue > obj.ContentDrivenSize.(direction)
                        obj.ContentDrivenSize.(direction) = childValue;
                    end

                case 'childRemoved'
                    if isempty(childLayout) ||...
                            (childLayout.(direction)(end) >= obj.ContentDrivenSize.(direction))
                        obj.ContentDrivenSize.(direction) = obj.computeContentDrivenLengthForSingleDirection(direction, child);
                    end
                    
                case 'recalculateAll'
                    obj.ContentDrivenSize.(direction) = obj.computeContentDrivenLengthForSingleDirection(direction);
            end
        end
        
        function [rows,cols] = computeImplicitGridSize(obj, action, varargin)
            % Compute the number of rows and columns needed to fit all
            % components with as few implicit rows and columns as possible
            %
            % If there are more components than can fit in the user-set
            % dimensions (RowHeight, ColumnWidth), then this calculates
            % and returns the number of implicit rows and columns needed to
            % contain all the components
            %
            % Else, if the user dimensions are greater than or equal to the
            % content-driven dimensions, this returns the user dimensions
            %
            % Example: If a user has a 2x2 grid with 4 components, and then
            % adds a 5th component, the grid will add an implicit row to
            % make space for the new component. This will be reflected in
            % the RowHeight, and the grid dimensions will become 3x2.
            %
            % If the user removes the added component, the implicit row
            % will be removed, returning the grid to its last explicitly
            % set dimensions of 2x2.
            %
            % Input: 
            %        -action -> reason for recomputing the content driven size; 
            %                   allowed: 'childAdded','childRemoved','recalculateAll'
            %        -varagin -> child component that is either being added or removed
            
            
            obj.updateContentDrivenLengthForSingleDirection('Row', action, varargin{:});
            obj.updateContentDrivenLengthForSingleDirection('Column', action, varargin{:});
            
            rows = obj.computeImplicitSizeForSingleDirectionFromLength(obj.LastUserSetSizes.RowHeight, obj.ContentDrivenSize.Row);
            cols = obj.computeImplicitSizeForSingleDirectionFromLength(obj.LastUserSetSizes.ColumnWidth, obj.ContentDrivenSize.Column);

        end

        
        % Update Implicit Grid Size
        function updateImplicitGridSize(obj, action, varargin)
            % Input: 
            %        -action -> reason for recomputing the content driven size; 
            %                   allowed: 'childAdded','childRemoved','recalculateAll'
            %        -varagin -> child component that is either being added or removed
            
            % Update the number of rows and columns in grid
            [rows,cols] = obj.computeImplicitGridSize(action, varargin{:});
            
            % Update object and send to the view
            if ~isequal(obj.PrivateRowHeight, rows)
                obj.PrivateRowHeight = rows;
                obj.markPropertiesDirty({'RowHeight'});
            end
            
            if ~isequal(obj.PrivateColumnWidth,cols)
                obj.PrivateColumnWidth = cols;
                obj.markPropertiesDirty({'ColumnWidth'});
            end
        end

        function updateLastCell(obj, action, varargin)
            % Update the LastCell property to reflect the latest
            %
            % INPUTS:
            % - action: either 'childAdded', 'childRemoved',
            %           'LayoutChanged' or 'ColumnWidthChanged'
            % - varargin: child being added or removed, if applicable
            
            if nargin == 3
                child = varargin{1};
                childLayout = child.Layout;
            end
            
            switch action
                case 'childAdded'            
                    if childLayout.Row(end) > obj.LastCell(1) || ...
                       ( childLayout.Row(end) >= obj.LastCell(1) && childLayout.Column(end) > obj.LastCell(2) )
                        % Component added after last cell
                        obj.LastCell = [childLayout.Row(end), childLayout.Column(end)];
                    end
                    
                case 'childRemoved'
                    obj.LastCell = obj.computeLastCell(child);
                
                case 'LayoutChanged'
                    obj.LastCell = obj.computeLastCell();
                
                case 'ColumnWidthChanged'
                    if ~obj.isLastCellValid()                        
                        obj.LastCell = obj.computeLastCell();
                    end
            end
        end
        
        function lastCell = computeLastCell(obj, varargin)
            % Computes the last occupied cell 'from scratch'.
            % 'Last' when scanning the cells left to right, then top to
            % bottom
            %
            % INPUT: 
            % - varargin: child being removed, if applicable
            
            gridChildren = obj.Children;
            
            if nargin == 2
                childRemoved = varargin{1};
                
                % Remove the child being deleted from the list.
                % The child being removed should not be accounted for
                % in determining the last cell
                index = gridChildren == childRemoved;
                gridChildren(index) = [];
            end
            
            lastCell = obj.getInitialLastCell();
            
            numChildren = length(gridChildren);
            
            for k = 1:numChildren
                child = gridChildren(k);
                
                childLayout = child.Layout;
                        
                % Row and Column can be a 1x2 array (spanning).
                % Use the last element to capture the furthest cell.
                row = childLayout.Row(end);
                col = childLayout.Column(end);

                if(row > lastCell(1))
                    % this cell is on a row further down, mark it as
                    % the new lastCell
                    lastCell = [row, col];
                end

                if(row == lastCell(1))
                    % the component is on the last row, check the column
                    if(col > lastCell(2))
                        lastCell(2) = col;
                    end
                end

            end
        end
        
        function initLastCell = getInitialLastCell(obj)
            
            % Initialize in the last cell of the row above the first row.
            % 
            % This is so we don't need to special case the scenarios where
            % - the grid has no children
            % - the grid is of size 0x0
            initLastCell = [0, length(obj.ColumnWidth)];
        end
        
        function isValid = isLastCellValid(obj)
            % Whether obj.LastCell is valid
            %
            % LastCell becomes invalid when ColumnWidth changes and there
            % are no children because LastCell is inialized based on
            % ColumnWidth
            % 
            % When there are children in the grid, LastCell is always valid
            % because of the implicit rows and columns
            
            initLastCell = obj.getInitialLastCell();
            
            if obj.LastCell(1) == initLastCell(1) && obj.LastCell(2) ~= initLastCell(2)
                % Grid is empty. Initial last cell value is invalid now
                % because the number of columns has changed
                isValid = false;
            else 
                isValid = true;
            end
        end
        
        function processChildAdded(obj, newChild)
            % A new child has been added to the container, figure out where
            % to place it in the grid
            %
            
            currentLayoutOptions = newChild.Layout;
            if obj.validateChildsLayoutOptions(currentLayoutOptions)
                return;
            end
            
            nextCell = obj.findNextAvailableCell();
            
            constraints = matlab.ui.layout.GridLayoutOptions;
            constraints.Row = nextCell(1);
            constraints.Column = nextCell(2);
            
            % Set Layout so the property can be marked dirty and
            % sent to the view
            newChild.setLayoutFromLayoutContainer(constraints);
            
        end
    end
    
    
    methods(Access = 'private', Static)
        
        function validateConstraintsPvPairs(pvPairs)
            % Validate that pvPairs has an even number of elements and that
            % all the property names are constraints that can be set for
            % the grid
            
            assert(mod(length(pvPairs), 2) == 0);
            
            propertyNames = pvPairs(1:2:end);
            
            % Get the cosntraints that are valid for the grid
            constraintsNames = properties('matlab.ui.layout.GridLayoutOptions');
            
            % Check that all property names are valid constraints
            isPropertyNameValid = ismember(propertyNames, constraintsNames);
            assert(all(isPropertyNameValid));            
        end
        
        function output = validateRowColumnSize(input)
            % Validate that input is a valid RowHeight or ColumnWidth value
            %
            % Input must be a row cell array where each element is either:
            % - a positive number (pixel value)
            % - 'fit' 
            % - 'Nx' where N is any positive number (weight)
            %
            % Strings are accepted and converted into char array
            % If the input is a column, it is accepted and converted into a
            % row
            %
            % e.g. 
            % input = {10, 'fit'};
            % input = {'fit', '1x', '3x'};
            %
            % Arrays are accepted and converted into cells
            % 
            % e.g. 
            % ["1x", "2x", "fit"] is converted to {'1x', '2x', 'fit'}
            % [100, 200, 300] is converted to {100, 200, 300}
            
            
            if(iscell(input) && isempty(input))
                % Treat {} separately because it doesn't pass the test of
                % 'vector' in validateattribute
                output = input;
                return
            end
            
            % Convert string array to cell of chars
            if(isstring(input))
                input = cellstr(input);
            end
            
            % Convert numeric array to cell of numbers
            if(isnumeric(input))
                input = num2cell(input);
            end
            
            % Verify that it is a vector 
            validateattributes(input, ...
                {'cell'}, ...
                {'vector'});
            
            % Reshape to row
            input = matlab.ui.control.internal.model.PropertyHandling.getOrientedVectorArray(input, 'horizontal');
            
            % Validate each element of the cell
            output = input;
            isValid = true;
            
            for k = 1:length(input)
                el = input{k};
                
                isPositiveNumber = matlab.ui.control.internal.model.PropertyHandling.isPositiveNumber(el);
                
                if ~isPositiveNumber
                    
                    el = convertStringsToChars(el);
                    
                    isRowCharArray = ischar(el) && isrow(el); 
                    if ~isRowCharArray
                        isValid = false;
                        break;
                    end
                    
                    isFit = strcmp(el, 'fit');
                    if isFit                        
                        isValid = true;                        
                    else
                        % check for '1x', '2x', ...
                        doesEndWithX = length(el) >= 2 && strcmp(el(end), 'x');
                        weight = str2double(el(1:end-1));
                        isWeightValid = matlab.ui.control.internal.model.PropertyHandling.isPositiveNumber(weight);
                        
                        isOfFormNx = doesEndWithX && isWeightValid;                        
                        if ~isOfFormNx
                            isValid = false;
                            break;
                        end
                        
                        % use weight directly to filter additional
                        % characters from the stored values:
                        % '+', leading zeros, trailing zeros
                        % Stored values will look like result of str2num
                        el = char(weight + "x");
                    end
                end
                
                % this element is valid, store it
                output{k} = el;                
            end
            
            if ~isValid                
                error('Invalid row/column size');                    
            end                  
        end
        
        
    end
    
end