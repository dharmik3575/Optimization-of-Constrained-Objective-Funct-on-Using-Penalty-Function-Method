%% constrv Checks for constraint voilation
% constrv(XVector,n)
%
% n == 1 is used to check if point is in feasible reagion
% n == 3 is used to return constraint violation
% n ~= (1 or 3) is used to return the penalty for 
% function violation by bracket Penatly

function out = constrv(x,n)
global problem
if problem == 1
    x1 = x(1);
    x2 = x(2);
    
    c(1) = min((x1-5)^2 + (x2-5)^2 - 100,0);
    c(2) = min(82.81 - ((x1-6)^2 + (x2-5)^2),0);
    c(3) = min((x1-13),0);
    c(4) = min((20-x1),0);
    c(5) = min(x2,0);
    c(6) = min(4-x2,0);
    
elseif problem == 2
    x1 = x(1);
    x2 = x(2);
    
    c(1) = min(-(x1^2 - x2 + 1),0);
    c(2) = min(-(1 - x1 + (x2-4)^2),0);
    c(3) = min(10-x1,0);
    c(4) = min(10-x2,0);
    c(5) = min(x1,0);
    c(6) = min(x2,0);
    
elseif problem == 3
    x1 = x(1);
    x2 = x(2);
    x3 = x(3);
    x4 = x(4);
    x5 = x(5);
    x6 = x(6);
    x7 = x(7);
    x8 = x(8);
    
    c(1) = min(1-0.0025*(x4+x6),0);
    c(2) = min(1-0.025*(-x4+x5+x7),0);
    c(3) = min(1-0.01*(-x6+x8),0);
    c(4) = min(-(100*x1 + (-x6*x1) + (833.33252*x4) + (-83333.333)),0);
    c(5) = min(-(x2*x4 + (-x2*x7) + (-1250*x4) + (1250*x5)),0);
    c(6) = min(-(x3*x5 + (-x3*x8) + (-2500*x5) + 1250000),0);
    c(7) = min(10000-x1,0);
    c(8) = min(10000-x2,0);
    c(9) = min(10000-x3,0);
    c(10) = min(x1-100,0);
    c(11) = min(x2-1000,0);
    c(12) = min(x3-1000,0);
    c(13) = min(1000-x4,0);
    c(14) = min(1000-x5,0);
    c(15) = min(1000-x6,0);
    c(16) = min(1000-x7,0);
    c(17) = min(1000-x8,0);
    c(18) = min(x4-10,0);
    c(19) = min(x5-10,0);
    c(20) = min(x6-10,0);
    c(21) = min(x7-10,0);
    c(22) = min(x8-10,0);
    
end

if n == 1
    out = logical(sum(c));
elseif n == 3
    out = sum(c);
else
    out = sum(c.^2);
end

end