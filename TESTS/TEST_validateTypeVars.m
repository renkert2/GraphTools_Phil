clear all
close all


T = Type('Val_Char',"xt");

% list of symbolic variables
T_var_all = symvar(T.Val_Sym).';

% list of options to define a symbolic expresssion
% note u is a vector of length = the number of unique variables. This
% should capture cases where a type is defined without state arguments
% ex: think if T = "u1*u2*u3*u4"...

T_var(ismember(T_var,[sym('xh');sym('xt')])) = []; % remove head and tail states from the list of options
uOptions = sym('u',[length(T_var) 1]); % total number of allowable inputs

% check to see if T_var is a subset of variable options
if sum(ismember(T_var,uOptions)) == length(T_var)
    disp('valid type')
else
    error('invalid type')
end


