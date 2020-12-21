classdef GraphVertex < matlab.mixin.Heterogeneous & matlab.mixin.Copyable
    %GRAPHVERTEX Super Class for all vertex types (i.e. State,
    %External)
    %   The matlab.mixin.Heterogeneous type will allow us to create heterogeneous object
    %   arrays of all the child classes.  
    properties
        Description string = string.empty()
        Type = 1 % Vertex Type: 1 - Energy Flow, 2 - State Flow
        Capacitance (:,1) Type_Capacitance = Type_Capacitance.empty();
        Coefficient (:,1) {mustBeNonnegative} = 0
        Initial = 0
    end
    
    methods
        function obj = GraphVertex(varargin)
            if nargin > 1
                obj = my_inputparser(obj,varargin{:});
            end
        end
    end
    
%     methods (Sealed)
%         function x = eq(A,B)
%             x = zeros(size(B));
%             for i = 1:length(B)
%                 try
%                     x(i) = A == B(i);
%                 catch
%                     x(i) = logical(0);
%                 end
%             end           
%         end             
%     end

end

