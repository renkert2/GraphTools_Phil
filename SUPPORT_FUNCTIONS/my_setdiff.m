function [ X ] = my_setdiff( A,B )

% more computationally efficient version of setdiff
check = false(1, max(max(A), max(B)));
check(A) = true;
check(B) = false;
X = A(check(A));


end

