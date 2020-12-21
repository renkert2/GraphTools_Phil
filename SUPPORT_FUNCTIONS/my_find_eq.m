function [x] = my_find_eq(A,B,offset)

% try 
    x = find(A == B)+offset;
    if isempty(x)
        x = 0;
    end
% catch
%     x = 0;
% end

end

