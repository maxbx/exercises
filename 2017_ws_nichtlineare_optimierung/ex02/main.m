function main

% function definition
f = @(x) function_name(x);
x1min = -5; x1max = 5; x2min = -5; x2max = 5;
plotfunction(f, x1min, x1max, x2min, x2max);

% random starting point
x01 = rand(1)*(x1max-x1min)+x1min;
x02 = rand(1)*(x2max-x2min)+x2min;
cpoint = [x01; x02];

figure(1);
scatter(cpoint(1), cpoint(2), 'k', 'filled');

[cval, cgrad] = f(cpoint);
while true
	% steepest descent with an random angle
	tau = pi/2*rand(1)-pi/4;
	dir = -[cos(tau) -sin(tau); sin(tau) cos(tau)]*cgrad;
	maxa = 1; c1 = 0.0001; c2 = 0.9;

	[npoint, nval, ngrad] = linesearch(f, cpoint, dir, maxa, cval, cgrad, c1, c2);

	% plot step of linesearch
	figure(1);
	N = 50;
	x = linspace(cpoint(1), npoint(1), N);
	y = linspace(cpoint(2), npoint(2), N);
	plot(x, y, 'b', 'Linewidth', 2); % line from cpoint to npoint
	scatter(npoint(1), npoint(2), 'k', 'filled'); % npoint

	% plot values along the line
	figure(2);
	clf;
	N = 50;
	alpha = linspace(0, maxa, N);
	x = cpoint+alpha.*dir;
	y = zeros(1, N);
	for ii = 1:length(x)
		y(ii) = f(cpoint+alpha(ii)*dir);
	end
	plot(alpha, y); %real values
	hold on;
	a = 1/(3*maxa^3)*(6*cval+3*maxa*cgrad'*dir-6*nval+3*maxa*ngrad'*dir);
	b = 1/maxa^2*(-3*cval-2*maxa*cgrad'*dir+3*nval-maxa*ngrad'*dir);
	c = cgrad'*dir;
	d = cval;
	plot(alpha, a*alpha.^3+b*alpha.^2+c*alpha+d); % interpolation
	%plot(alpha, c1*cgrad'*dir*alpha+y(1)); %Armijo
	%plot(alpha, c2*cgrad'*dir*alpha+y(1)); %Curvature
	calpha = sqrt(sum((npoint-cpoint).^2))./sqrt(sum(dir.^2)); % alpha for npoint
	scatter(calpha, f(cpoint+calpha*dir), 'filled');
	drawnow;
	pause

	if abs(cval-nval) < 10^-2 % if progress is little end
		return;
	end

	cpoint = npoint;
	cval = nval;
	cgrad = ngrad;
end
end


function [val, grad] = function_name(x)
%
% gives back the value, the gradient and the hessian of a function
%

% Rosenbrock
%val = 100*(x(2, :)-x(1, :).^2).^2-(1-x(1, :)).^2;
%grad = [
%	-400*(x(2, :).*x(1, :)-x(1, :).^3)+2-2*x(1, :);
%	200*(x(2, :)-x(1, :).^2)
%];

% another function
val = sin(3*x(1, :))+sin(x(2, :));
grad = [3*cos(3*x(1, :)); cos(x(2, :))];

end


