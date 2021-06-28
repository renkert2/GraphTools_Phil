classdef System < SystemElement
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
        Components (:,1) SystemElement
        ConnectP (:,1) cell
        
        Model Model = Model.empty()
    end
    
    methods
        function obj = System(Name, Comps, arg3)
            if nargin == 1
                obj.Name = Name;
            end
            if nargin == 2
                obj.Name = Name;
                obj.Components = Comps;
            end
            if nargin == 3
                obj.Name = Name;
                obj.Components = Comps;
                
                if isa(arg3, 'cell')
                    % Create system from components and connections
                    obj.ConnectP = arg3;
                    obj.init_super();
                    obj=SystemCombine(obj);
                elseif isa(arg3, 'Model')
                    % Assign predefined model
                    obj.Model = arg3;
                end
            end
        end
        
        function obj = SystemCombine(obj)
            obj = Combine(obj.Components, obj.ConnectP, obj);
        end
        
        function model = createModel(obj, varargin)
            % Creates GraphModel from obj.Graph and stores
            % it in obj.Model.  Pass additional arguments 
            % to the GraphModel constructor with varargin
            
            model = GraphModel(obj, varargin{:});
            obj.Model = model;
        end
        
        function model = get.Model(obj)
            % Get method creates a GraphModel of obj.Graph
            % if one doesn't already exist
            
            if isempty(obj.Model)
                createModel(obj);
            end
            model = obj.Model; 
        end
        
        function C = getComponents(obj, Type)
            i = arrayfun(@(x) isa(x,Type), obj.Components);
            C = obj.Components(i);
        end
    end
end

