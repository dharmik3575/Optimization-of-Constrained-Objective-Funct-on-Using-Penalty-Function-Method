%% Marquart's Method
% Marquart(XVector,RValue,MaximumItterations)
%
% Calculates the extrema by Marquart's Method

function res = Marquart(x,R,M)

global feval problem
fprintf("\n");
disp("R = "+num2str(R))

gd = grad(x,R);                                 % Gradient of the function at initial Guess
k = 0;                                          % Itteration Counter
lam = min(10^length(char(string(R)))*1e5,1e10); % Lamda value depends on R
fx = func(x,R,1);                               % Evaluating Function Value at Initial Guess
check = true;                                   % Check to stop while loo
c2 = 10;                                        % Check for direction if are dependent or independent
c1 = 6;                                         % Check For Gradient if very low
red = 0;                                        % Check if restart is required

while check
    if sum(isnan(x))~=0
        error("x vector contains NaN. Try Again with different Initial Guess");
    end
    res(k+1,1) = k;                             % res is a matrix that stores data
    res(k+1,2) = fx;
    res(k+1,3) = func(x,R,2);
    res(k+1,4) = constrv(x,3);
    res(k+1,5) = feval;
    for ij = 1:length(x)                        % Storing the vector at itteration 'k'
        res(k+1,ij+5) = x(ij);
    end
    
    H = hess(x,R);
    s = (H + lam*eye(size(x,1)));
    if (((norm(gd)<(1/(10^c1)))||(k>=M)))       % Checking if termination conditions are satisfied
        check = false;
        fprintf("Gradient Low ===> TERMINATING\n");
        
    elseif (abs(det(s))<(1e-6))||(rcond(s)<1e-10)                 % Restart if determinant of hessian is very low
        for ii = 1:size(x,1)
            rng('shuffle');
            x(ii) = x(ii)+randi(ceil(abs(max(x)))+randi(5,1,1),1,1);
        end
        red=0;
        fprintf("Restart Performed\n");
        fx = func(x,R,1);
        
    else
        s = -1*(s\gd);                               % New Direction
        s = s*UniD(x,s,R);                           % UniD performes the unidirectional search
        
        x = x + s         ;                          % New value
        if(k~=0)
            c2 = abs(s1.'*s);
            %c2 = norm(x - x1)/norm(x1);
        end
        if((c2<1/(R*5))||(c2<1e-6))&&(k>1)           % Check if Directions are Independent
            check = false;
            fprintf("Directions not Independent ===> TERMINATING\n");
        else
            k = k+1;
            fx1 = func(x,R,1);
            s1 = s;
            if fx1<fx                                % Check to increase or decrease Lambda Value
                lam = max(lam/10,0.001);
            else
                lam = lam*10;
            end
            
            if((abs(fx1-fx)<1e-6)&&(k>1))            % Doing a restart if there is not much improvement in funciton value
                red = red+1;
                if red==5
                    for ii = 1:size(x,1)
                        rng('shuffle');
                        x(ii) = x(ii)*randi(10,1,1);
                    end
                    red=0;
                end
                fprintf("Restart Performed\n");
                fx1 = func(x,R,1);
            end
            
            fx = fx1;
            gd = grad(x,R);
        end
    end
end

% Storing Data in a Table for better visualization of the data
var = {};
var{1} = 'Itteration';
var{2} = 'Penalty_Func_val';
var{3} = 'Function_Val';
var{4} = 'Contraint_Violation';
var{5} = 'Function_Evaluations';
for ii=6:size(res,2)
    var{ii} = char("x" + string(ii-4));
end
res = array2table(res,'VariableNames',var);

end

%% Hessian Matrix
% hess(x_vector,RValue)
%
% Hessian Matrix is calculated numerically

function out = hess(x,R)
a = size(x,1);
out = zeros(a);
for ii = 1:a
    for jj = 1:a
        out(ii,jj) = pardef(x,[ii jj].',R);
    end
end

out = 00.5*(out+out.');                 % This is to make the hessian matrix symmetric
end

%% Gradient of the function
% grad(x_vector,RValue)
%
% Gradient is calculated numerically

function out = grad(x,R)
out = zeros(size(x,1),1);
for ii=1:size(x,1)
    out(ii) = pardef(x,ii,R);
end
end

%% Partial differentiation
% pardef(x_vector,xj,Rvalue)
%
% This performs the partial differentiation
% numerically which is required to calculate
% The hessian and gradient

function out = pardef(x,m,R)
if size(m,1) == 1
    h = (1e-4)*m;
    xm = x;
    xn = x;
    xm(m) = x(m) + h;
    xn(m) = x(m) - h;
    out = (func(xm,R,1) - func(xn,R,1))/(2*h);

elseif size(m,1) == 2
    h1 = x*0;
    h1(m(1)) = 1e-2*m(1);
    h2 = x*0;
    h2(m(2)) = 1e-2*m(2);
    out = (func(x+h1+h2,R,1) - (func(x+h1,R,1) + func(x+h2,R,1)) + func(x,R,1))/(h1(m(1))*h2(m(2)));
end
end
