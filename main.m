%% Constraint Optimization Using Penalty Function Method And Marquardt Method
%% Reading the Problem Number

clear;clc;
global problem feval
fp = fopen("input.txt",'r');
problem = str2double(fgetl(fp));    % 'problem' is the problem number
fclose(fp);
%% Processing to solution

fprintf("Initial Guess ===>\n");
x = randvec(problem);                % Initial Guess
disp(x)

tic;                                % Calculating Time
Rmax = 7;                           % Maximum Itteration For Penalty Function
M = 100;                            % Maximum Itteration For Marquardt Method
result = PenaltyFunc(x,7,100);      % Penalty Function method
a = toc;
%% Printing and Ploting Results
fprintf("\n-----------------------------------------------\n")
fprintf("\tResult by the program\n")
fprintf("-----------------------------------------------\n")
fprintf("\n\nProblem Number = %d\n",problem);
fprintf("Time Taken = %fs\n",a);
res = table2array(result{1,2});
fprintf("\n\nInitial Guess ===>\n");
disp(res(1,6:end).');
fprintf("\nInitial Value of R = %d\n",result{1,1});
fprintf("Penalty Function Value = %f\n",res(1,2));
fprintf("Initial Value of Function = %f\n",res(1,3));
fprintf("Constraint Violation = %f\n",res(1,4));
figure(1)
res = result{1,2};
res = res(:,1:5);
stackedplot(res,'XVariable','Itteration');
title("Convergence Plot (R = "+num2str(result{1,1})+")")

res = table2array(result{end,2});
fprintf("\n\nFinal Vector ====>\n");
disp(res(end,6:end).');
fprintf("\nFinal Value of R = %d\n",result{end,1});
fprintf("Penalty Function Value = %f\n",res(end,2));
fprintf("Final Value of Function = %f\n",res(end,3));
fprintf("Constraint Violation = %f\n",res(end,4));
fprintf("Total Function Evaluations = %d\n",feval);
figure(2)
res = result{end,2};
res = res(:,1:5);
stackedplot(res,'XVariable','Itteration');
title("Convergence Plot (R = "+num2str(result{end,1})+")")

figure(3)
res = result{1,2};
for ii = 2:size(result,1)
    res = [res; result{ii,2}]; 
end
res = res(:,2:5);
res.Itterations = [1:height(res)].';
stackedplot(res,'XVariable','Itterations');
title("Combined Convergence Plot for all values of R");
fprintf("\n-----------------------------------------------\n")
fprintf("-----------------------------------------------\n")
%% Contour Plot

figure(4)
if problem ~=3
a = result{1,2};
xv(1,1:2) = table2array(a(1,6:7));
for ii = 1:size(result,1)
   a = result{ii,2};
   xv(ii+1,1:2) = table2array(a(end,6:7));
end
a = ceil(max(max(abs(xv))));
[X, Y] = meshgrid(-a:0.01:a,-a:0.01:a);
if problem == 1
    F = (X-10).^3 + (Y-20).^3;
    opt = [14.095 0.84296];
elseif problem == 2
    F = (((sin(2*pi*X)).^3).*sin(2*pi*Y))./((X.^3).*(X+Y));
    opt = [1.227 4.225];
end
contour(X,Y,F);

hold on
plot(xv(:,1),xv(:,2),'LineWidth',2)
plot(opt(1),opt(2),'*')
legend({'Contour','by procedure','Actual Optimum'})
xlabel("x1 --->")
ylabel("x2 --->")
title("Contour Plot for Problem "+num2str(problem))
hold off
end

%% Results by fmincon
% This is used to compare results with fmincon function.
% 
% To change value of 'R' change the value in 'funcfmin.m' all other perimeters 
% are not required to be changed. R == 1e6 unless changed.
% 
% Uncomment below code to use fmincon


objective = @funcfmin;
fprintf("\n\n-----------------------------------------------\n")
fprintf("\tResults by fmincon ===>\n")
fprintf("-----------------------------------------------\n")
fprintf("\nInitial Function Value = %f\n",func(x,0,2));

[x] = fmincon(objective,x,[],[],[],[]);
fprintf("\nFinal vector ===>\n")
disp(x)
fprintf("\nFinal Function value ===>\n")
disp(func(x,0,2))
fprintf("Constraint Violation ===>\n")
disp(constrv(x,3))
fprintf("\n-----------------------------------------------\n")
fprintf("-----------------------------------------------\n")


%% Writing Data to file
% output is saved to 'OUTPUT.mat' in form of cell data type.

save('OUTPUT.mat','result');
%% Function to Select Random Vector
% This will be the initial guess. The initial guess is choosen infeasible that 
% is checked by function 'constrv'

function out = randvec(x)
if x == 1
    x = 2;
elseif x==2
    x = 2;
elseif x == 3
    x = 8;
end
for ii = 1:50
    rng("shuffle");
    for jj = 1:x
    out = randn(x,1)*(ii+7);
    end
    if constrv(out,1)
        break
    end
end
end
