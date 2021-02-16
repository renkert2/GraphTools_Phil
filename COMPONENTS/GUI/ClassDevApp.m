classdef ClassDevApp < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure           matlab.ui.Figure
        TabGroup           matlab.ui.container.TabGroup
        ComponentTab       matlab.ui.container.Tab
        GraphTab               matlab.ui.container.Tab
        PropertiesPanel    matlab.ui.container.Panel
        GraphGroup    matlab.ui.container.TabGroup
        VertexTab       matlab.ui.container.Tab
        EdgeTab               matlab.ui.container.Tab
        
        GridLayout         matlab.ui.container.GridLayout
        AddPropertyButton  matlab.ui.control.Button
        BuildClassButton  matlab.ui.control.Button
        DeletePropertyButton (:,1)  matlab.ui.control.Button = matlab.ui.control.Button.empty()
        PropertyLabel      matlab.ui.control.Label
        DefaultValueLabel  matlab.ui.control.Label
        PropertyCommentLabel matlab.ui.control.Label
        PropertyFields  (:,1) MultiField = MultiField.empty()
    end

    % Callbacks that handle component events
    methods (Access = private)

        % Button pushed function: AddPropertyButton
        function AddPropertyButtonPushed(app, event)
            idx = app.AddPropertyButton.Layout.Row;
            n = 3;
            app.PropertyFields(end+1) = MultiField(app.GridLayout,n);
            app.GridLayout.RowHeight  = [app.GridLayout.RowHeight, 20];
            app.AddPropertyButton.Layout.Row = idx+1;     
            
            app.DeletePropertyButton(end+1) = uibutton(app.GridLayout,'Text','Delete Property');
            app.DeletePropertyButton(end).Layout.Row = idx;
            app.DeletePropertyButton(end).Layout.Column = 1;
            app.DeletePropertyButton(end).ButtonPushedFcn = createCallbackFcn(app, @DeletePropertyButtonPushed, true);
        end
        
        % Button pushed function: AddPropertyButton
        function DeletePropertyButtonPushed(app, event)
            delIdx = app.DeletePropertyButton == event.Source;
            delete(app.PropertyFields(delIdx).Fields); app.PropertyFields(delIdx) = [];
            delete(app.DeletePropertyButton(delIdx)); app.DeletePropertyButton(delIdx) = [];
            for i = find(delIdx):numel(app.PropertyFields)
                for j = 1:numel(app.PropertyFields(i).Fields)
                    app.PropertyFields(i).Fields(j).Layout.Row = app.PropertyFields(i).Fields(j).Layout.Row-1;
                end
                app.DeletePropertyButton(i).Layout.Row = app.DeletePropertyButton(i).Layout.Row-1;
            end
            app.AddPropertyButton.Layout.Row = app.AddPropertyButton.Layout.Row-1;
            app.GridLayout.RowHeight = app.GridLayout.RowHeight(1:end-1);
        end
        
        function BuildClassButtonPushed(app, event)
            Name = 'Test';
            AllFields = [app.PropertyFields(:).Fields];
            Props = {AllFields(1,:).Value};
            Val = {AllFields(2,:).Value};
            Com = {AllFields(3,:).Value};
            CusCom = 'None';
            
            BuildCompClass(Name,Props,Val,Com,CusCom)
            

        end
        
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create UIFigure and hide until all components are created
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.Position = [100 100 689 480];
            app.UIFigure.Name = 'UI Figure';

            % Create TabGroup
            app.TabGroup = uitabgroup(app.UIFigure);
            app.TabGroup.Position = [253 1 437 480];

            % Create ComponentTab
            app.ComponentTab = uitab(app.TabGroup);
            app.ComponentTab.Title = 'Component';
            
            % Create Tab2
            app.GraphTab = uitab(app.TabGroup);
            app.GraphTab.Title = 'GraphModel';

            % Create TabGroup
            app.GraphGroup = uitabgroup(app.GraphTab);
            app.GraphGroup.Position = [1 100 436 355];

            % Create ComponentTab
            app.VertexTab = uitab(app.GraphGroup);
            app.VertexTab.Title = 'Vertices';
            
            % Create Tab2
            app.EdgeTab = uitab(app.GraphGroup);
            app.EdgeTab.Title = 'Edges';
            
            
            
            % Create PropertiesPanel
            app.PropertiesPanel = uipanel(app.ComponentTab);
            app.PropertiesPanel.Title = 'Properties';
            app.PropertiesPanel.Position = [1 100 436 355];

            % Create GridLayout
            app.GridLayout = uigridlayout(app.PropertiesPanel);
            app.GridLayout.ColumnWidth = {100, '1x', '1x','1x'};
            app.GridLayout.RowHeight = {20, 20};

            % Create AddPropertyButton
            app.AddPropertyButton = uibutton(app.GridLayout, 'push');
            app.AddPropertyButton.ButtonPushedFcn = createCallbackFcn(app, @AddPropertyButtonPushed, true);
            app.AddPropertyButton.Layout.Row = 2;
            app.AddPropertyButton.Layout.Column = 1;
            app.AddPropertyButton.Text = 'Add Property';

            % Create BuildClassButton
            app.BuildClassButton = uibutton(app.UIFigure, 'push');
            app.BuildClassButton.ButtonPushedFcn = createCallbackFcn(app, @BuildClassButtonPushed, true);
            app.BuildClassButton.Text = 'BuildClass';
            app.BuildClassButton.Position= [10 10 100 22];
        
            % Create PropertyLabel
            app.PropertyLabel = uilabel(app.GridLayout);
            app.PropertyLabel.HorizontalAlignment = 'center';
            app.PropertyLabel.Layout.Row = 1;
            app.PropertyLabel.Layout.Column = 2;
            app.PropertyLabel.Text = 'Property Name';

            % Create DefaultValueLabel
            app.DefaultValueLabel = uilabel(app.GridLayout);
            app.DefaultValueLabel.HorizontalAlignment = 'center';
            app.DefaultValueLabel.Layout.Row = 1;
            app.DefaultValueLabel.Layout.Column = 3;
            app.DefaultValueLabel.Text = 'Default Value';

            % Create DefaultValueLabel
            app.PropertyCommentLabel = uilabel(app.GridLayout);
            app.PropertyCommentLabel.HorizontalAlignment = 'center';
            app.PropertyCommentLabel.Layout.Row = 1;
            app.PropertyCommentLabel.Layout.Column = 4;
            app.PropertyCommentLabel.Text = 'Comment';
            

            % Show the figure after all components are created
            app.UIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = ClassDevApp

            % Create UIFigure and components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.UIFigure)

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.UIFigure)
        end
    end
end