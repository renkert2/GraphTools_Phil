classdef LinearModel < Model
    % LinearModel is a propety of the model class and has linear state
    % space model matrices. The model is of the form:
    %
    % System Description
    % x_dot = A*x + B*u + E*d + f0
    % y     = C*x + D*u + H*d + g0
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Contributors: Christopher T. Aksland and Phil Renkert
    % Association: University of Illionis at Urbana-Champaign
    % Contact: aksland2@illinois.edu and renkert2@illinois.edu
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % potential improvements:
    % f0 and g0 should be moved into the linear model class...
    % Finish CalcState and CalcOutput
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    properties      
        A_sym  (:,:) {mustBeNumericOrSym}
        B_sym  (:,:) {mustBeNumericOrSym}
        E_sym  (:,:) {mustBeNumericOrSym}
        
        C_sym  (:,:) {mustBeNumericOrSym}
        D_sym  (:,:) {mustBeNumericOrSym}
        H_sym  (:,:) {mustBeNumericOrSym}
        
        f0_sym (:,1) {mustBeNumericOrSym}
        g0_sym (:,1) {mustBeNumericOrSym}
    end
        
    properties (Access = private)   
        A_func  function_handle {mustBeScalarOrEmpty}
        B_func  function_handle {mustBeScalarOrEmpty}
        E_func  function_handle {mustBeScalarOrEmpty}
        
        C_func  function_handle {mustBeScalarOrEmpty}
        D_func  function_handle {mustBeScalarOrEmpty}
        H_func  function_handle {mustBeScalarOrEmpty}
        
        f0_func function_handle {mustBeScalarOrEmpty}
        g0_func function_handle {mustBeScalarOrEmpty}
    end
    
    properties (Constant, Access = private)
        N_mats = 6
        syms_props = ["A_sym", "B_sym", "E_sym", "C_sym", "D_sym", "H_sym", "f0_sym", "g0_sym"];
        func_props = ["A_func", "B_func", "E_func", "C_func", "D_func", "H_func", "f0_func", "g0_func"];
    end

    methods
        function obj = LinearModel(varargin) % constructor stores the linear matrix information
            if nargin
                switch nargin
                    case 6
                        obj.A_sym = varargin{1};
                        obj.B_sym = varargin{2};
                        obj.E_sym = varargin{3};
                        obj.C_sym = varargin{4};
                        obj.D_sym = varargin{5};
                        obj.H_sym = varargin{6};
                        
                        obj.f0 = zeros(size(obj.A_sym, 2),1);
                        obj.g0 = zeros(size(obj.C_sym, 1),1);
                    case 8
                        obj.A_sym = varargin{1};
                        obj.B_sym = varargin{2};
                        obj.E_sym = varargin{3};
                        obj.C_sym = varargin{4};
                        obj.D_sym = varargin{5};
                        obj.H_sym = varargin{6};
                        obj.f0 = varargin{7};
                        obj.g0 = varargin{8};
                end
                obj.init();
            end
        end
        
        function init(obj)
            obj.Nx = size(obj.A_sym, 2);
            obj.Nu = size(obj.B_sym, 2);
            obj.Nd = size(obj.E_sym, 2);
            obj.Ny = size(obj.C_sym, 1);
            
            if isempty(obj.SymVars)
                setSymVars(obj);
            end
            
            setLinearCalcFuncs(obj)
            setFGSym(obj);
            setCalcFuncs(obj)
        end
        
        function setLinearCalcFuncs(obj)
            for i = 1:obj.N_mats
                sym = obj.(obj.syms_props(i));
                if ~isempty(sym)
                    func = genMatlabFunctions(obj,sym);
                    obj.(obj.func_props(i)) = func;
                end
            end
        end
        
        function setFGSym(obj)
            x = obj.SymVars.x;
            u = obj.SymVars.u;
            d = obj.SymVars.d;
            
            obj.f_sym = sum([obj.A_sym*x, obj.B_sym*u, obj.E_sym*d, obj.f0_sym],2); % using horzcat and sum allows us to ignore empty values
            obj.g_sym = sum([obj.C_sym*x, obj.D_sym*u, obj.H_sym*d, obj.g0_sym],2);
        end
            
        function [A,B,E,C,D,H] = CalcMatrices(obj,x,u,d)
            param_lengths = [obj.Nx, obj.Nu, obj.Nd];

            vars = {x,u,d};
            
            for i = 1:numel(vars)
                assert(size(vars{i},1) >= param_lengths(i), "Argument %d requires %d entries", i, param_lengths(i));
            end
            
            mats = cell(obj.N_mats,1);
            for i = 1:obj.N_mats
                func = obj.(obj.func_props(i));
                if ~isempty(func)
                    mats{i} = obj.CalcX(func, vars);
                end
            end
            [A,B,E,C,D,H] = mats{:};
        end
    end

    methods (Static)
        function [A_d, B_d, C_d, D_d] = Discretize(A,B,C,D,dT,opts)
            % Discretizes a Linear State Space Model, see details on
            % Wikipedia page
            arguments
                A
                B
                C
                D
                dT (1,1) double
                opts.AbsTol double = 1e-8
                opts.RelTol double = 1e-8
            end

            A_d = expm(A*dT);

            fun = @(dt) expm(dt*A);
            Aint  = integral(fun,0,dT,'AbsTol',1e-8,'RelTol',1e-8,'ArrayValued', true);
            B_d = Aint*B;
            
            if ~isempty(C)
                C_d = C;
            end
            if ~isempty(D)
                D_d = D;
            end
        end
    end
end

