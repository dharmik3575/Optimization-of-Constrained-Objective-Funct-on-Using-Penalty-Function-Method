function out = funcfmin(x)
R = 1e6;
C = 1;
global problem
if problem == 1
    out = (x(1)-10)^3 + (x(2)-20)^3;
    %out = x(1)^2 + 2*x(2)^2;
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
