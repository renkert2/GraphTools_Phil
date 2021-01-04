function Yi = qminterp1(X,Y,xi)
% QMINTERP1 1-dimensional fast interpolation (for 1D vectors)
%
% Usage:
%   yi = qminterp1(X,Y,xi)  - Similar usage to interp1
%
% Usage restrictions
%   X must be monotonic and evenly spaced
%   X and Y must be VECTORS
%   xi is a scalar, not a vector
%   Only linear interpolation is used   
%   Presently, no extrapolation is performed
%
%
% Error checking
%   WARNING: Little error checking is performed on the X vector. If it is
%   not monotonic and evenly spaced, erroneous results will be
%   returned.
%
% Author: T.L. McKinley - 23 May '07
%

% find vector increment and lower bounds
delx = X(2)-X(1);
xref = (xi - X(1))/delx;
ixlow = floor(xref)+1;

if ( ixlow > 0 && ixlow < length(X) )
%
%  interpolate using natural coordinates
%
    xnat = (xi - X(ixlow))/delx;
    Yi = Y(ixlow) + (Y(ixlow+1) - Y(ixlow)) * xnat;
else
    Yi = NaN;
end
