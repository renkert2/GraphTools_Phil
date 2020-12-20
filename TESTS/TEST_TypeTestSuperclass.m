%% Capacitance Type
% t_cap = Type_Capacitance('a^2') % - Returns invalid type definition message

% syms a
% t_cap = Type_Capacitance(a^2) % - Returns invalid type definition message

t_cap = Type_Capacitance('x^2'); % - success

syms x
t_cap = Type_Capacitance(x^2); % - success

t_cap.Val_Str = "x^3";

t_cap.Val_Sym = x^5;

t_cap.calcVal(2)
t_cap.calcJac(3,3)
%% PowerFlow Type

t_pf = Type_PowerFlow('xt*xh*u3')

t_pf.Val_Str = ('xt*u7')

syms xt xh

t_pf.Val_Sym = xt*xh



%% Calc Power Flows:
t_pf = Type_PowerFlow('xt*xh*u3')
xt_ = 1;
xh_ = 2;
u_ = [1 2 3];
t_pf.calcVal(xt_, xh_, u_)

t_pf = Type_PowerFlow('xt*xh')
xt_ = 1;
xh_ = 2;
t_pf.calcVal(xt_, xh_)
