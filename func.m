%% Objective Function
% func(XVector,RValue,C)
%
% C == 1 for penaly function
% C ~= 1 for function value

function out = func(x,R,C)
global problem feval

feval = feval + 1;

if problem == 1
    out = (x(1)-10)^3 + (x(2)-20)^3;
    if C == 1
        out = out + R*constrv(x,2);
    end
    
elseif problem == 2
    out = (((sin(2*pi*x(1)))^3)*(sin(2*pi*x(2))))/(x(1)^3*(x(1)+x(2)));
    if C == 1
        out = -out + R*constrv(x,2);
    end
    
elseif problem == 3
    out = x(1) + x(2) + x(3);
    if C == 1
        out = out + R*constrv(x,2);
    end
end
end
