classdef paramFit < handle
    % Wraps Curve Fitting Toolbox fits used for data-based surrogate models
    % that use continuous functions to estimate dependent compParams from
    % other compParams
    properties
        Inputs (:,1) compParam
        Outputs (:,1) compParam
        
        BoundaryWarning logical = true
    end
    
    properties (SetAccess = private)
        Data (:,1) struct % Struct array (obj.N_outs, 1), each with Inputs and Outputs field for use with model creation
        Boundary Boundary
        Models (:,1) cell % Cell array of fit objects or function handles
        Functions (:,1) cell % Cell array of function handles, wrap Models with additional functionality like boundary checking
    end
    
    properties (SetAccess = private)
        N_ins double
        N_outs double
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
            
            for i = 1:obj.N_outs
                model = makeFit(obj, obj.Data(i).Inputs, obj.Data(i).Outputs, fit_type{i}, fit_opts{i});
                obj.Models{i} = model;
                obj.Functions{i} = makeModelWrapper(obj,model);
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
        
        function outs = calcParams(obj, varargin)
            outs = zeros(obj.N_outs,1);
            for i = 1:obj.N_outs
                f = obj.Functions{i};
                outs(i,1) =  f(varargin{:});
            end
        end
        
        function plot(obj, varargin)
            if nargin > 1
                outs = calcParams(obj, varargin{:});
            end
            
            for i = 1:obj.N_outs
                figure(i)
                plot(obj.Models{i}, obj.Data(i).Inputs, obj.Data(i).Outputs);
                
                title('paramFit Plot')
                [olb,oub] = bounds(obj.Data(i).Outputs);

                switch obj.N_ins
                    case 1
                        pfun = @plot;
                        ylim([olb oub]);
                        if ~isempty(obj.Inputs)
                            xlabel(latex(obj.Inputs(1)),'Interpreter','latex');
                        end
                        if ~isempty(obj.Outputs)
                            ylabel(latex(obj.Outputs(i)),'Interpreter','latex');
                        end
                    case 2
                        pfun = @plot3;
                        zlim([olb oub]);
                        if ~isempty(obj.Inputs)
                            xlabel(latex(obj.Inputs(1)),'Interpreter','latex');
                            ylabel(latex(obj.Inputs(2)),'Interpreter','latex');
                        end
                        if ~isempty(obj.Outputs)
                            zlabel(latex(obj.Outputs(i)),'Interpreter','latex');
                        end
                end
                
                if nargin > 1
                    hold on
                    pfun(varargin{:}, outs(i),'.r','MarkerSize',20)
                    hold off
                end
            end
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

