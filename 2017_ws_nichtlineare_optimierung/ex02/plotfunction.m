function plotfunction(f, x1min, x1max, x2min, x2max)
%
% plots a function
%

steps = 100;
[x1, x2] = meshgrid(linspace(x1min, x1max, steps), linspace(x2min, x2max, steps));

val = zeros(steps);
for n1 = 1:steps
	for n2 = 1:steps
		val(n1, n2) = f([x1(n1, n2); x2(n1, n2)]);
	end
end

% contourp plot
figure(1);
clf;
%contour(x1, x2, val, [0.1, 0.5, 2.5, 12.5, 62.5, 312.5]);
contour(x1, x2, val, linspace(min(val(:)), max(val(:)), 20));
hold on;

axis([x1min, x1max, x2min, x2max, -10, 1000]);
caxis([min(val(:)), max(val(:))]);

end
