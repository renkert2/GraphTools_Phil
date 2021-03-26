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
        
        f0 (:,1) {mustBeNumericOrSym}
        g0 (:,1) {mustBeNumericOrSym}
    end
        
    properties (Access = private)   
        A_func  function_handle {mustBeScalarOrEmpty}
        B_func  function_handle {mustBeScalarOrEmpty}
        E_func  function_handle {mustBeScalarOrEmpty}
        
        C_func  function_handle {mustBeScalarOrEmpty}
        D_func  function_handle {mustBeScalarOrEmpty}
        H_func  function_handle {mustBeScalarOrEmpty}
    end
    
    properties (Constant, Access = private)
        N_mats = 6
        syms_props = ["A_sym", "B_sym", "E_sym", "C_sym", "D_sym", "H_sym"];
        func_props = ["A_func", "B_func", "E_func", "C_func", "D_func", "H_func"];
    end

    methods
        function obj = LinearModel(varargin) % constructor stores the linear matrix information
            if nargin == 6
                obj.A_sym = varargin{1};
                obj.B_sym = varargin{2};
                obj.E_sym = varargin{3};
                obj.C_sym = varargin{4};
                obj.D_sym = varargin{5};
                obj.H_sym = varargin{6};
                
                obj.Nx = size(obj.A_sym, 2);
                obj.Nu = size(obj.B_sym, 2);
                obj.Nd = size(obj.E_sym, 2);
                obj.Ny = size(obj.C_sym, 1);
                
                obj.f0 = zeros(obj.Nx,1);
                obj.g0 = zeros(obj.Ny,1);
                
                obj.init();
            end
        end
        
        function init(obj)
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
            
            if isempty(u)
                u = 0;
            end
            if isempty(d)
                d = 0;
            end
            
            obj.f_sym = obj.A_sym*x + obj.B_sym*u + obj.E_sym*d + obj.f0;
            obj.g_sym = obj.C_sym*x + obj.D_sym*u + obj.H_sym*d + obj.g0;
        end
            
        function [A,B,E,C,D,H] = CalcMatrices(obj,x,u,d,params)
            param_lengths = [obj.Nx, obj.Nu, obj.Nd, obj.SymParams.N];
            
            if isempty(obj.SymParams)
                vars = {x,u,d};
            else
                if nargin == 4
                    vars = {x,u,d,obj.SymParams.Vals};
                elseif nargin == 5
                    vars = {x,u,d,params};
                end
            end
            
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
end

