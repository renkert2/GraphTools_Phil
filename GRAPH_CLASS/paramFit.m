classdef paramFit < handle
    % Wraps Curve Fitting Toolbox fits used for data-based surrogate models
    % that use continuous functions to estimate dependent compParams from
    % other compParams
    properties
        Inputs (:,1) compParam
        Outputs (:,1) compParam
        
        FitTypes (:,1) cell
        FitOpts (:,1) cell % cell of structs
        
        BoundaryWarning logical = true
    end
    
    properties (SetAccess = private)
        Data (:,1) struct % Struct array (obj.N_outs, 1), each with Inputs and Outputs field for use with model creation
        Boundary Boundary
        Models (:,1) cell % Cell array of fit objects or function handles
        Functions (:,1) cell % Cell array of function handles, wrap Models with additional functionality like boundary checking
    end
    
    properties (SetAccess = private)
        N_ins double = []
        N_outs double = []
    end
    
    methods
        function obj = paramFit(ins, outs)
            if nargin == 2
                if isa(ins, 'compParam') && isa(outs, 'compParam')
                    obj.Inputs = ins;
                    obj.Outputs = outs;
                elseif isnumeric(ins) && isnumeric(outs)
                    mustBeInteger(ins);
                    mustBeInteger(outs);
                    obj.N_ins = ins;
                    obj.N_outs = outs;
                else
                    error("invalid constructor arguments.  Must be the input and output compParams arrays or the number of inputs and outputs")
                end
            end
        end
        
        function set.Inputs(obj, ins)
            obj.Inputs = ins;
            obj.N_ins = numel(ins);
        end
        
        function set.Outputs(obj, outs)
            obj.Outputs = outs;
            obj.N_outs = numel(outs);
        end
        
        function set.FitTypes(obj, types)
            if obj.N_outs
                assert(numel(types) == obj.N_outs, "Each output must be associated with a fit")
            end
            obj.FitTypes = types;
        end
        
        function set.FitOpts(obj, opts)
            if obj.N_outs
                assert(numel(opts) == obj.N_outs, "Each output must be associated with a fit")
            end
            obj.FitOpts = opts;
        end
        
        function setData(obj, input_data, output_data)
            N_points = size(input_data,1);
            assert(N_points == size(output_data,1), 'Input and Output Data must have equal number of rows');
            
            for i = 1:obj.N_outs
                [obj.Data(i,1).Inputs, obj.Data(i,1).Outputs] = prepareData(obj, input_data, output_data(:,i));
            end
        end
        
        function setBoundary(obj)
            obj.Boundary = Boundary(obj.Data(1).Inputs);
        end
        
        function setModels(obj, fit_type, fit_opts)
            arguments
                obj
            end
            arguments (Repeating)
                fit_type
                fit_opts
            end
            if nargin > 1
                obj.FitTypes = fit_type;
                obj.FitOpts = fit_opts;
            end
            
            for i = 1:obj.N_outs
                model = makeFit(obj, obj.Data(i).Inputs, obj.Data(i).Outputs, obj.FitTypes{i}, obj.FitOpts{i});
                obj.Models{i} = model;
                obj.Functions{i} = makeModelWrapper(obj,model);
            end
            
            if ~isempty(obj.Outputs)
                setOutputDependency(obj);
            end
            
            function fh = makeModelWrapper(obj,model)
                fh = @ModelWrapper;
                
                function out = ModelWrapper(varargin)
                    if obj.BoundaryWarning
                        if all(isnumeric([varargin{:}]))
                            in_bounds = obj.Boundary.isInBoundary(varargin{:});
                            if any(~in_bounds)
                                out_i = find(~in_bounds);
                                out_str = num2str(out_i, '%d, ');
                                warning("Points %s Outside Boundary", out_str)
                            end
                        end
                    end
                    out = model(varargin{:});
                    % Can add post process functionality like LB < out < UB later
                end
            end
        end
        
        function setOutputDependency(obj)
            for i = 1:obj.N_outs
                f = obj.Functions{i};
                setDependency(obj.Outputs(i), f, obj.Inputs);
            end
        end
        
        function setSpan(obj, span)
           assert(isnumeric(span) && span >= 0 && span <= 1, "Span must be between 0 and 1");
           for i = 1:obj.N_outs
               ft = obj.FitTypes{i};
               if checkL(ft)
                   obj.FitOpts{i}.Span = span;
               end
           end
           setModels(obj);
           
            function i = checkL(ft)
               t = string(type(ft));
               i = (t == "lowess" | t == "loess");
            end
        end
        
        function outs = calcParams(obj, varargin)
            outs = zeros(obj.N_outs,1);
            for i = 1:obj.N_outs
                f = obj.Functions{i};
                outs(i,1) =  f(varargin{:});
            end
        end
        
        function F = plot(obj, opts)
            arguments
                obj
                opts.Outputs = 1:obj.N_outs
                opts.PlotMarker logical = false
            end
            
            if ~isempty(obj.Inputs) && ~isempty(obj.Outputs)
                ins = {obj.Inputs.Value};
                outs = calcParams(obj, ins{:});
            else
                ins = {};
                outs = [];
            end
            
            out_plts = opts.Outputs;
            F = matlab.ui.Figure.empty();            
            for i = 1:numel(out_plts)
                f = figure(i);
                F(i) = f;
                
                out_I = out_plts(i);
                plot(obj.Models{out_I}, obj.Data(out_I).Inputs, obj.Data(out_I).Outputs);
                colormap()
                
                [olb,oub] = bounds(obj.Data(out_I).Outputs);

                switch obj.N_ins
                    case 1
                        pfun = @plot;
                        ylim([olb oub]);
                        if ~isempty(obj.Inputs)
                            xlabel(latex(obj.Inputs(1)),'Interpreter','latex');
                        end
                        if ~isempty(obj.Outputs)
                            ylabel(latex(obj.Outputs(out_I)),'Interpreter','latex');
                        end
                    case 2
                        pfun = @plot3;
                        zlim([olb oub]);
                        if ~isempty(obj.Inputs)
                            xlabel(latex(obj.Inputs(1)),'Interpreter','latex');
                            ylabel(latex(obj.Inputs(2)),'Interpreter','latex');
                        end
                        if ~isempty(obj.Outputs)
                            zlabel(latex(obj.Outputs(out_I)),'Interpreter','latex');
                        end
                end
                
                if ~isempty(outs) && opts.PlotMarker
                    hold on
                    pfun(ins{:}, outs(out_I),'or', 'MarkerSize',10, 'LineWidth', 2, 'MarkerEdgeColor', 'w', 'MarkerFaceColor','r')
                    hold off
                end
            end
        end
        
        function F = plotErrors(obj, opts)
            arguments
                obj
                opts.Outputs = 1:obj.N_outs
            end
            
            out_plts = opts.Outputs;
            F = matlab.ui.Figure.empty();
            
            x = linspace(0,1,256/2)';
            R = [1 0 0];
            W = [1 1 1];
            B = [0 0 1];
            map = [B+x.*(W - B); W+x.*(R - W)];
            
            for i = 1:numel(out_plts)
                f = figure(i);
                F(i) = f;
                
                out_I = out_plts(i);
                
                input_data = obj.Data(out_I).Inputs;
                output_data = obj.Data(out_I).Outputs;
                model = obj.Models{out_I};
                model_data = model(input_data); % Evaluate the model at the design points
                error = (model_data - output_data)./output_data;
                
                
                scatter(input_data(:,1),input_data(:,2),[],error,'filled',...
                  'MarkerEdgeColor', 0.8*[1 1 1],...
                  'LineWidth',1)
                
                c = colorbar;
                c.Label.String = 'Rel. Error';
                c.Label.Interpreter = 'latex';
                c.Label.FontSize = 12;
                c.FontSize = 12;
                colormap(map)
                
                padding = 0.1;
                
                xminmax = minmax(input_data(:,1)');
                xrange = range(xminmax);
                xlim(xminmax + padding.*xrange.*[-1,1])
                
                yminmax = minmax(input_data(:,2)');
                yrange = range(yminmax);
                ylim(yminmax + padding.*yrange.*[-1,1])
              
                me = max(abs(error));
                caxis([-me, me]);
                
                if ~isempty(obj.Inputs)
                    xlabel(latex(obj.Inputs(1)),'Interpreter','latex');
                    ylabel(latex(obj.Inputs(2)),'Interpreter','latex');
                end

            end
        end
        
        function p = plotBoundary(obj, varargin)
            plot(obj.Boundary, varargin{:});
            switch obj.N_ins
                case 1
                    xlabel(latex(obj.Inputs(1)),'Interpreter','latex');
                case 2
                    xlabel(latex(obj.Inputs(1)),'Interpreter','latex');
                    ylabel(latex(obj.Inputs(2)),'Interpreter','latex');
            end
            p = gca;
        end
        
        function cftool(obj, output_number)
            dat = obj.Data(output_number);
            switch obj.N_ins
                case 1
                    cftool(dat.Inputs(:,1),dat.Outputs);
                case 2
                    cftool(dat.Inputs(:,1),dat.Inputs(:,2),dat.Outputs);
            end
        end
    end
    
    methods (Access = protected)
        function [inputs_out, outputs_out] = prepareData(obj, inputs_in, outputs_in)
            % Override in subclasses for custom data prep
            switch obj.N_ins
                case 1
                    [inputs_out, outputs_out] = prepareCurveData(inputs_in, outputs_in);
                case 2
                    [inputs_out(:,1), inputs_out(:,2), outputs_out] = prepareSurfaceData(inputs_in(:,1),inputs_in(:,2), outputs_in);
                otherwise
                    inputs_out = inputs_in;
                    outputs_out = outputs_in;
            end
        end
        
        function [model] = makeFit(obj, input_data, output_data, ft, fo)
            if obj.N_ins == 1 || obj.N_ins == 2
                model = fit(input_data, output_data, ft, fo);
            else
                error("paramFit.Models not compatible with more than two input variables")
            end
        end
    end
end

