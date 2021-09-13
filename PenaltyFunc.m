%% Penalty Function method
% PenaltyFunc(Initual_Guess,Maximum_Itr_R,Max_Itr_Marquardt)
%
% Uses Marquardts Method to get extrema value

function res = PenaltyFunc(xvec,n,M)

global feval
feval = 0;                      % Feval stores function evaluations

c = 10;
R = c;

% 'res' is cell Data Type It stores the results
% First Column contains 'R' Values
% Second Column Contains table of results by 
% Marquardt Method for each itteration
res{1,1} = R;
res{1,2} = Marquart(xvec,R,M);

xvec = res{end,2};
disp("Itterations = "+num2str(height(xvec)-1));
xvec = table2array(xvec(end,6:end)).';
fprintf("x vector ===>\n")
disp(xvec);
for ii=2:n
    if sum(isnan(xvec))~=0
        error("x vector contains NaN. Try Again with different Initial Guess");      % Error is there is NaN in the x Vector
    end
    R = c*R;
    res{ii,1} = R;
    res{ii,2} = Marquart(xvec,R,M);
    a = res{ii,2};
    b = res{ii-1,2};
    xvec = res{end,2};
    disp("Itterations = "+num2str(height(xvec)-1));
    xvec = table2array(xvec(end,6:end)).';
    fprintf("x vector ===>\n")
    disp(xvec);
    
    if abs(a.Penalty_Func_val(end) - b.Penalty_Func_val(end))<1e-4
        fprintf("Penalty Function Method Terminated\n");
        break;
    end
end
end
