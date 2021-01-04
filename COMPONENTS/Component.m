classdef Component < matlab.mixin.Heterogeneous & handle
    %COMPONENT Super class to be inherited by Components (i.e.
    %Battery, Motor, Heat Exchanger, etc...
    %   Detailed explanation goes here
    
    properties %(SetAccess = protected)        
        % Block Name
        Name string = "Component"
        graph Graph = Graph.empty()
    end
    
    methods
        function obj = Component(varargin)
            if nargin == 1
                if isstruct(varargin{1})
                    for i = 1:numel(varargin{1})
                        try
                            obj.(varargin{1}(i).Name) = varargin{1}(i).Value;
                        catch
                            % eventually update this to indicate that no
                            % property for this class exists.
                        end
                    end
                else
                    error('Components must be defined using a structure of Name/Value fields or Name-Value pairs')
                end
            elseif nargin > 1
                my_inputparser(obj,varargin{:}); % input parser component models
            end
            obj.init(); % I don't know why we need this and can't just call ConstructGraph - CTA

        end
        
        function set.Name(obj, name)
            obj.Name = string(name);
        end
        
        function init(obj)
            obj.ConstructGraph();
            obj.DefineChildren();
        end
    end
    
    methods (Access = protected)
        function ConstructGraph(obj)
            g = DefineGraph(obj);
            obj.graph = g;
        end
        
        function DefineChildren(obj)
            try
                for i = 1:numel(obj.graph.Inputs)
                    obj.graph.Inputs(i).Parent = obj;
                end
            end
        end
    end
    
    methods (Abstract, Access = protected)
        DefineGraph(p)       
    end
end

