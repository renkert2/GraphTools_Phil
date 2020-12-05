classdef Component < matlab.mixin.Heterogeneous
    %GRAPHCOMPONENT_SUPER Super class to be inherited by Components (i.e.
    %Battery, Motor, Heat Exchanger, etc...
    %   Detailed explanation goes here
    
    properties (SetAccess = protected)
        GraphModel GraphModel
    end
    
    methods
        function obj = GraphComponent(p)
            if nargin == 1
                props = fieldnames(p);
                for iprop = 1:length(props)
                    try
                        obj.(props{iprop}) = p.(props{iprop});
                    catch
                        warning('No Property %s in component', props{iprop});
                    end
                end
                
                obj.init();
            end
        end
        
        function init(obj)
            obj.ConstructGraphModel();
        end
    end
    
    methods (Access = protected)
        function ConstructGraphModel(obj)
            g = DefineGraphModel(obj);
            g.init();
            obj.GraphModel = g;
        end
        
        function g = DefineGraphModel(p)
            g = GraphModel(); % Function to be defined by child classes
        end
    end
    
    
end

