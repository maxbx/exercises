function main
%
% plots the Rosenbrock function, the linear and the quadratic model of a random point
%

f = @(x1, x2) rosenbrock(x1, x2);
x1min = -1; x1max = 2; x2min = -1/2; x2max = 3;

plotfunction(f, x1min, x1max, x2min, x2max);

x1 = rand(1)*(x1max-x1min)+x1min;
x2 = rand(1)*(x2max-x2min)+x2min;

linear_model(f, [x1, x2], x1min, x1max, x2min, x2max);
quadratic_model(f, [x1, x2], x1min, x1max, x2min, x2max);

end
