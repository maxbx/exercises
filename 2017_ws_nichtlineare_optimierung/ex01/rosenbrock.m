function [val, grad, hessian] = rosenbrock(x1, x2)
%
% gives back the value, the gradient and the hessian of the rosenbrock function
%

val = 100*(x2-x1.^2).^2-(1-x1).^2;

grad = [
	-400*(x2.*x1-x1.^3)+2-2*x1,
	200*(x2-x1.^2)
];

hessian = [
	-400*x2+1200*x1.^2-2, -400*x1;
	-400*x1, 200*ones(size(x1))
];

end
