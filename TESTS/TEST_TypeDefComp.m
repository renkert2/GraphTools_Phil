clear all
close all

N = 20;
method = 2;

%% Method of defining few power flows with coeffs
for j = 1:10
if method == 1
    tic
    T(1) = Type_PowerFlow("xt*xh*u1");
    T(2) = Type_PowerFlow("xt^2*u1");
    T(3) = Type_PowerFlow("xt-xh");
    
    for i = 1:N
        idx = max(mod(i,4),1);
        E1(i) = GraphEdge_Internal('PowerFlow',T(idx),'Coefficient',idx);
    end
    
    PTypeAll = vertcat(E1.PowerFlow); % list of all capacitance types
    numPType = arrayfun(@(x) length(x.PowerFlow),E1); % find number of capacitance types per vertex
    [P_coeff,PType] = MakeCoeffMatrix(E1,PTypeAll,numPType);
    
    t(j) = toc;
    
    %% Method of defining Unique power flows with coeffs
    
elseif method == 2
    tic
    for i = 1:N
        idx = max(mod(i,4),1);
        switch idx
            case 1
                E2(i) = GraphEdge_Internal('PowerFlow',Type_PowerFlow("1*xt*xh*u1"));
            case 2
                E2(i) = GraphEdge_Internal('PowerFlow',Type_PowerFlow("2*xt^2*u1"));
            case 3
                E2(i) = GraphEdge_Internal('PowerFlow',Type_PowerFlow("3*(xt-xh)"));
        end
        
    end
    t(j) = toc;
end
end

% the first call seems to be the slowest, so lets just omit it...
% however, in general the first call is the most important because we will
% have to call this functions a first time.
mean(t(2:end))
std(t(2:end))

