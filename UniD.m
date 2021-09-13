%% Unidirection Search
% UniD(Xvector,Direction,RValue)
%
% Tries to perform unidirection search in given direction
% If Newton Raphson do not converge but region is given by bounding phase
% method then the mid-point of the region is returned. In this case Region
% will not be large so returned value will be close to optimum.

function out = UniD(x,s,R)
[a,b] = bounding_phase(min(norm(s),randi(10)),x,s,R,100000);
[out, check] = NR(a,b,x,s,R);
if check ~=0
    [a,b] = bounding_phase(out,x,s,R,1000);
    [out, check] = NR(a,b,x,s,R);
    if check ~=0
        [a,b] = bounding_phase(out,x,s,R,1);
        out = (a+b)/2;
    end
end
end

%% Bounding Phase Method
% bounding_phase(initial_guess,XVector,Direction,RValue,MaxRegionSize)
%
% Bounding Phase return the region in which the optimum point is possible
% The maximum size or the region is given by MaxRegionSize which is
% increased in case region is not found even after 3700 itterations.

function [a, b] = bounding_phase(ini,x,s,R,checker)
delta = 0.01;

% This for loop is to ensure that search is started in the right direction
% in case the Intial Guess do not lead to any direction it is incremented
% or decremented by 5.

for ii = 1:5
a3 = ini - delta;
a4 = ini + delta;
f3 = func(x+a3*s,R,1);
f0 = func(x+ini*s,R,1);
f4 = func(x+a4*s,R,1);

if f4 > f0
    if f0 > f3
        delta = (-1) * delta;
        break;
    else
        ini = ini+5;
    end
elseif f3 > f0
    if f0 > f4
        delta = 1*delta;
        break;
    else
        ini = ini-5;
    end
end
end

a1 = ini;
k = 0;
itr = 0;
check = true;
Nmax = 5000;
m = 2^k * delta;
while check
    a2 = a1 + m;
    f1 = func(x+a1*s,R,1);
    f2 = func(x+a2*s,R,1);
    if f1<f2
        check = false;
    elseif itr>Nmax
        error("Error in Unidirection Search. Coudn't bound by Bounding Phase Method. Try with any other point");
    else
        ini = a1;
        a1 = a2;
        k = k + 1;
        itr = itr +1;
        
        if m > checker                      % If Region becomes large 'k' is made zero
            if itr>=Nmax*0.7
                checker = checker*10;       % Increasing Region size if point is not bounded
            end
            k = 0;
        end
        
        m = 2^k * delta;
    end
end

if ini < a2
    a = ini; b = a2;
else
    a = a2; b = ini;
end
end

%% Newton Raphson
% NR(a,b,XVector,Direction,RValue)
%
% Newton Raphson return the optimum value in the region given by the
% bounding phase method. In case the method do not converge it return the
% midpoint of the region.
% checker == 0 ===> Newton Raphson Converged
% checker == (1 or 2) ===> Newton Raphson Did not converge

function [x_next, checker] = NR(a,b,x,s,R)
x1 = a;

ii = 0;
Nmax = 500;
d = 0.001;
e = 1e-4;

check = 1;
red = 0;

checker = 0;
while(check == 1)
    f_x = Diff(x1,1,d,R,x,s);
    x_next = x1 - (f_x/Diff(x1,2,d,R,x,s));
    if((x_next>b)||(x_next<a))
        if red == 0
            red = 1;
            x1 = b;
        elseif red == 1
            red = 2;
            x1 = (a+b)/2;
        else
            x_next = (a+b)/2;
            checker = 1;
            break;
        end
        f_x = Diff(x1,1,d,R,x,s);
    end
    if abs(f_x)<e
        check = 0;
    elseif ii >= Nmax
         x_next = (a+b)/2;
         checker = 2;
         break;
    else
        x_next = x1 - (f_x/Diff(x1,2,d,R,x,s));
        ii = ii+1;
        x1 = x_next;
    end
end
end

%% Derivative
% Diff(alpha,n,deltaValue,RValue,XVector,Direction)
%
% Calculating Derivative by central difference method
% n == 1 ===> First Derivative
% n == 2 ===> Second Derivative

function out = Diff(a,n,del,R,x,s)
if n ==1
    out = (func(x+(a+del)*s,R,1) - func(x+(a-del)*s,R,1))/(2*del);
elseif n==2
    out = (func(x+(a+del)*s,R,1) - (2*func(x+a*s,R,1)) + func(x+(a-del)*s,R,1))/(del*del);
end
end