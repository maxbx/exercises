function quadratic_model(f, x0, x1min, x1max, x2min, x2max)
%
% calculates and plots the quadratic model of f in x0
%

[val, grad, hessian] = f(x0(1), x0(2));

steps = 50;
[x1, x2] = meshgrid(linspace(x1min, x1max, steps), linspace(x2min, x2max, steps));

figure(1);
contour(x1, x2, val+grad(1)*(x1-x0(1))+grad(2)*(x2-x0(2))+1/2*(hessian(1, 1)*(x1-x0(1)).^2+hessian(1, 2)*(x1-x0(1)).*(x2-x0(2))+hessian(2, 1)*(x2-x0(2)).*(x1*x0(1))+hessian(2, 2)*(x2-x0(2)).^2));
scatter(x0(1), x0(2), 'k', 'filled');

figure(2);
mesh(x1, x2, val+grad(1)*(x1-x0(1))+grad(2)*(x2-x0(2))+1/2*(hessian(1, 1)*(x1-x0(1)).^2+hessian(1, 2)*(x1-x0(1)).*(x2-x0(2))+hessian(2, 1)*(x2-x0(2)).*(x1*x0(1))+hessian(2, 2)*(x2-x0(2)).^2));
scatter3(x0(1), x0(2), f(x0(1), x0(2)), 'k', 'filled');



end
