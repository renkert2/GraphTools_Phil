classdef Component < matlab.mixin.Heterogeneous & handle
    %COMPONENT Super class to be inherited by Components (i.e.
    %Battery, Motor, Heat Exchanger, etc...
    %   Detailed explanation goes here
    
    properties %(SetAccess = protected)
        Model GraphModel
    end
    
    methods
        function obj = Component(varargin)
            if nargin > 1
                my_inputparser(obj,varargin{:}); % input parser component models
                obj.init(); % I don't know why we need this and can't just call ConstructGraph - CTA
            end
        end
        
        function init(obj)
            obj.ConstructGraph();
        end
    end
    
    methods (Access = protected)
        function ConstructGraph(obj)
            g = DefineGraph(obj);
            g.init();
            obj.Model = g;
        end
        
        function g = DefineGraph(p)
            g = GraphModel(); % Function to be defined by child classes
        end
    end
    
    
end

