function [p_min, min, iter] = pattern_search(f, p0, step_size, tolerance, N_it)

% Pattern search method
% It finds the relative minimum of a two variables function with the
% Pattern search method.
% The inputs are:
% - f: matlab inline function of the function that you want minimize;
% - p0: initial point for the search algorithm
% - step_size: initial step used for the exploration phase
% - tolerance: it defines a stopping criterion on the minimum size of 
%              exploration phase's step_size
% - N_it: it defines a stopping criterion on the max iterations number
% The outpus are:
% - iter: the number of iterations done
% - min: the function's value in the relative minimum finded.
% - p_min: it's a vector for the minimum point of the objective function
%
% IMPORTANT: for a correct work, the function 'pattern_search' must be preceded by a contour plot of the
% objective function to minimize.

f1=f;                       % Objective function
p(1,:)=p0;                  % Initial point
x0=p(1,:);
f=f1(p(1,1),p(1,2));        % Objective function's value in the initial point
f_min=f;
p_min=x0;                   
step=step_size;             % Step size is 0.5e-6
N=N_it;                     % Number of max iteration
l=2;                        % Support index
x = linspace(-10, 10, 20);
y = linspace(-10, 10, 20);
[X,Y] = meshgrid(x,y);
contour(x,y,f1(X,Y));
hold on;

plot(p(1,1),p(1,2),'ob')    % Initial point plot

xlim([-10 10]) 
ylim([-10 10]) 
iter=0;                     % Number of iteration done
ex_iter=0;                  % Number of successive extrapolations failed
ex_trov=0;                  % It's 1 when extrapolation has successfully
toll=tolerance;             % Stopping criterion on step-size
n_var=2;                    % Number of objective function's variables 

while iter<N && step>toll    
    spost=zeros(n_var*2,n_var);     % Construction of the movements matrix (contains transfers to all directions)    
    for k=1:n_var
        spost((k*2)-1,k)=step;      
        spost(k*2,k)=-step;
    end    
    k=1;                    % It's the index to change the direction
    trovato=0;              % It's 1 when the exploration has successfully    
    while k<=n_var*2   
    % Exploration
    p(l,:)=p(l-1,:)+spost(k,:);     % New point
    f(l)=f1(p(l,1),p(l,2));         % Objective function's value in the new point
    if f(l)>=f_min
        plot(p(l,1),p(l,2),'xr')    % Plot a 'x' if the point is worse than the precedent point
        plot([p(l-1,1) p(l,1)],[p(l-1,2) p(l,2)], 'LineStyle', ':', 'LineWidth', 1, 'Color', 'r')
        xlim([-10 10]) 
        ylim([-10 10]) 
        pause(0.2);
        k=k+1              % Increments the index for the movements matrix
    else
        plot(p(l,1),p(l,2),'og')    % Plot a 'o' if the point is better than the precedent point   
        plot([p(l-1,1) p(l,1)],[p(l-1,2) p(l,2)], 'LineStyle', '-', 'LineWidth', 1, 'Color', 'g')
        xlim([-10 10]) 
        ylim([-10 10]) 
        pause(0.2);
        f_min=f(l);         % Update the current function's minimum value 
        p_min=p(l,:);       % Update the current minimum point
        l=l+1;              % Increases the vector's index of the point finded 
        trovato=1;          % The exploration successfully
        if mod(k,2)==0      % If the exploration has found a new minimum point, it changes the exploration's direction 
            k=k+1
        else
            k=k+2
        end
    end
    end
    
    if ex_trov==1 || trovato==1 && ex_iter>0
        ex_iter=ex_iter-1;
    end 
    
    % It checks the successive extrapolation failed
    if ex_iter==1
        l=l-1;              % Reset
        x0=p(l-1,:);        % It's the new initial point for the extrapolation
        ex_iter=0;
        step=step/2;        % Decreases the exploration step-size
        f_min=f(l-1);
        ex_trov=0;
    else if trovato==0 && ex_iter==0 && ex_trov==0    % If, after a reset, failure occurs the extrapolation doesn't execute, but the step size is decreased
                step=step/2;             % Decreases the step-size         
                else     
                    % Extrapolation
                    ex_trov=0;
                    p(l,:)=p(l-1,:)+(p(l-1,:)-x0);
                    x0=p(l-1,:);
                    f(l)=f1(p(l,1),p(l,2));
                    if f(l)<f_min
                        f_min=f(l);
                        ex_trov=1;                        
                    end
                    plot(p(l,1),p(l,2),'sc')
                    plot([p(l-1,1) p(l,1)],[p(l-1,2) p(l,2)], 'LineStyle', '-', 'LineWidth', 2, 'Color', 'c')
                    xlim([-10 10]) 
                    ylim([-10 10]) 
                    pause(0.2);
                    l=l+1;
                    ex_iter=ex_iter+1;             
         end
    end
    iter=iter+1
end

disp('Il minimo della funzione si ha per: (var1,var2)')
disp(p_min)
disp('Il minimo della funzione ?: ')
min=f_min;
disp(f_min)