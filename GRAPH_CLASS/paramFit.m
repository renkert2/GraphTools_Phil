classdef paramFit < handle
    % Wraps Curve Fitting Toolbox fits used for data-based surrogate models
    % that use continuous functions to estimate dependent compParams from
    % other compParams
    properties
        Inputs (:,1) compParams
        Outputs (:,1) compParams
    end
    
    properties (SetAccess = private)
        Data (:,1) struct % Struct array (obj.N_outs, 1), each with Inputs and Outputs field for use with model creation
        Boundary Boundary
        Models (:,1) cell % Cell array of fit objects or function handles
    end
    
    properties (Access = private)
        N_ins double
        N_outs double
    end
    
    methods
        function obj = paramFit(ins, outs)
            obj.Inputs = ins;
            obj.Outputs = outs;
        end
        
        function set.Inputs(obj, ins)
            obj.Inputs = ins;
            obj.N_ins = numel(inputs);
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
            obj.Boundary = Boundary(obj.Data.Inputs);
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
                obj.Models{i} = makeFit(obj, obj.Data(i).Inputs, obj.Data(i).Ouptuts, fit_type{i}, fit_opts{i});
            end 
        end
                        
        function setOutputDependency(obj)
        end
        
        function outs = calcParams(obj, ins)
            outs = zeros(obj.N_outs,1);
            for i = 1:obj.N_outs
                f = obj.Models{i};
                outs(i,1) =  f(
        end
        
        function plot(obj)
        end
    end
    
    methods (Access = protected)
        function [inputs_out, outputs_out] = prepareData(obj, inputs_in, outputs_in)
            % Override in subclasses for custom data prep
            switch obj.N_ins
                case 1
                    [inputs_out, outputs_out] = prepareCurveData(inputs_in, outputs_in);
                case 2
                    [inputs_out, outputs_out] = prepareSurfaceData(inputs_in, outputs_in);
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

