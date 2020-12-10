clear all
close all


T = Type('Val_Char',"xt*xh*u1*u2*u3");

% list of symbolic variables
T_var = symvar(T.Val_Sym).';

% list of options to define a symbolic expresssion
% note u is a vector of length = the number of unique variables. This
% should capture cases where a type is defined without state arguments
% ex: think if T = "u1*u2*u3*u4"... 
varOptions = ['xt';'xh';sym('u',[length(T_var) 1])];

% check to see if T_var is a subset of variable options
if sum(ismember(T_var,varOptions)) == length(T_var)
    disp('valid type')
else
    error('invalid type')
end


% NOTE: varOptions could probably be reduced in size if we first detect whether
% states are included in the type definition. Then, we can exactly detect how many 
% inputs can be included

