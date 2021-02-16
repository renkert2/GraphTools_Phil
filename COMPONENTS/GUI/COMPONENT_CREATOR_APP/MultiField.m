classdef MultiField < handle 


    properties
        Fields (:,1) matlab.ui.control.EditField
%         Button (1,1) matlab.ui.control.Button
%         a = 0
    end


    
    methods
        function obj = MultiField(Parent,n)
            for i = 1:n
                obj.Fields(i) = uieditfield('Parent',Parent);
            end
%             obj.Button = uibutton('Parent',Parent,'Text','Delete Property','ButtonPushedFcn','@obj.DelPush');
        end
        
%         function DelPush(obj)
%             obj.a = obj.a + 1;
%         end
        
    end

end


