%%
% funCus = @(x)(20 + x(1)^2+x(2)^2-10*(cos(2*pi*x(1))+cos(2*pi*x(2))))
% lb = [-10,-10];
% ub = [10,10];
% A = [];
% b = [];
% Aeq = [];
% beq = [];
% x0 = [5,5];
% options = optimoptions('patternsearch','Display','iter','PlotFcn',@psplotbestf);
% nonlcon = [];
% x = patternsearch(funCus,x0,A,b,Aeq,beq,lb,ub, nonlcon, options)

%%
f  = inline('20 + x^2+y^2-10*(cos(2*pi*x)+cos(2*pi*y))', 'x', 'y');
x0 = [3.2, 5.8];
step_size = 1;
tolerance = 1e-8;
N_it = 1e5;
pattern_search(f, x0, step_size, tolerance, N_it)

%%
