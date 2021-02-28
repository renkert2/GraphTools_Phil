classdef System < handle
    % System is a superclass for any system composed of 
    % multiple components.
    
    % sys = System(Components, ConnectP)
    % - Create a System object from component array
    % Components and port connections defined in ConnectP.  See 
    % Component.Combine for more information
    
    % sys = System(Graph)
    % - Creates a System object from a precomposed Graph
    
    % Note: To save processing time, a GraphModel object is 
    % not created by default.  If the Model is empty
    % when System.Model is called, the get.Model method 
    % will construct a GraphModel with default property values.
    % If you need to pass additional arguments to the GraphModel 
    % constructor, use System.createModel(varargin) to pass 
    % varargin{:} to the GraphModel constructor
    
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
            if numel(varargin) > 0
                if numel(varargin) == 1
                    % Graph Passed to Constructor - Skip Initialization
                    obj.Graph = varargin{1};
                elseif numel(varargin) == 2
                    % Components and ConnectP Passed to Constructor - call init_super
                    obj.Components = varargin{1};
                    obj.ConnectP = varargin{2};
                    obj.init_super();
                end
            end
        end
               
        function init_super(obj)
            init(obj); % Optionally defined by subclasses
            createSysGraph(obj); % Calls Component.Combine to create system graph and resulting extrinsic props
        end
        
        function createSysGraph(obj)
            [obj.Graph, obj.extrinsicProps] = Combine(obj.Components, obj.ConnectP);
        end
        
        function model = createModel(obj, varargin)
            % Creates GraphModel from obj.Graph and stores
            % it in obj.Model.  Pass additional arguments 
            % to the GraphModel constructor with varargin
            
            obj.Model = GraphModel(obj.Graph, varargin{:});
            model = obj.Model;
        end
        
        function model = get.Model(obj)
            % Get method creates a GraphModel of obj.Graph
            % if one doesn't already exist
            
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

