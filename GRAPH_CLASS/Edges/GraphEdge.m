classdef GraphEdge < matlab.mixin.Heterogeneous & matlab.mixin.Copyable
    %GRAPHEdge Super Class for all edge types (i.e. Internal,
    %External)
    %   The matlab.mixin.Heterogeneous type will allow us to create heterogeneous object
    %   arrays of all the child classes.  
    
   


    methods
        function obj = GraphEdge(varargin)
            if nargin > 1
                obj = my_inputparser(obj,varargin{:});
            end
        end
    end

end

