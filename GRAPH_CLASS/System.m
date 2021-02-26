classdef System < handle
    % System ...
    % @Phil
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Contributors: Christopher T. Aksland and Phil Renkert
    % Association: University of Illionis at Urbana-Champaign
    % Contact: aksland2@illinois.edu and renkert2@illinois.edu
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % potential improvements:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    properties
        Name string = string.empty()
        
        Components (:,1) Component
        ConnectP (:,1) cell
        
        Graph Graph = Graph.empty()
        
        Model GraphModel = GraphModel.empty()
        Ports ComponentPort = ComponentPort.empty()
        
        extrinsicProps (:,1) extrinsicProp = extrinsicProp.empty()
    end
    
    methods
        function obj = System(varargin)
            if numel(varargin) == 0
                % Do Nothing
            elseif numel(varargin) == 1
                % Graph Passed to Constructor - Skip Initialization
                obj.Graph = varargin{1};
            elseif numel(varargin) == 2
                % Components and ConnectP Passed to Constructor - call init_super to 
                obj.Components = varargin{1};
                obj.ConnectP = varargin{2};
                obj.init_super();
            end
        end
               
        function init_super(obj)
            init(obj); % Optionally defined by subclasses
            createSysGraph(obj); % Calls Component.Combine to create system graph and resulting extrinsic props
        end
        
        function createSysGraph(obj)
            [obj.Graph, obj.extrinsicProps] = Combine(obj.Components, obj.ConnectP);
        end
        
        function model = createModel(obj)
            if isempty(obj.Model)
                obj.Model = GraphModel(obj.Graph);
            end
            model = obj.Model;
        end

        
        %% Methods to be overriden by subclasses
        function init(obj)
            % Placeholder
        end
    end
end

