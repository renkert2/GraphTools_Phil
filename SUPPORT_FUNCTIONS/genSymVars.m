function s = genSymVars(patt, n)
if n == 1
    s = sym(sprintf(patt,1));
elseif n>1
    s = sym(patt, [n 1]);
else
    s = sym.empty();
end
end
