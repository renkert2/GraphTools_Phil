% t_cap = Type_Capacitance('a^2') % - Returns invalid type definition message

% syms a
% t_cap = Type_Capacitance(a^2) % - Returns invalid type definition message

t_cap = Type_Capacitance('x^2'); % - success

syms x
t_cap = Type_Capacitance(x^2); % - success

t_cap.Val_Str = "x^3";

t_cap.Val_Sym = x^5;


