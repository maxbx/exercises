function plotfunction(f, x1min, x1max, x2min, x2max)
%
% plots a function
%

steps = 50;
[x1, x2] = meshgrid(linspace(x1min, x1max, steps), linspace(x2min, x2max, steps));


% contourp plot
figure(1);
clf;
contour(x1, x2, f(x1, x2), exp(-2:10)-min(f(x1,x2)(:)));
hold on;

axis([x1min, x1max, x2min, x2max, -10, 1000]);
caxis([0, 500]);


% mesh plot
figure(2);
clf;
mesh(x1, x2, f(x1, x2));
hold on;

axis([x1min, x1max, x2min, x2max, -10, 1000]);
caxis([0, 500]);
view(3);

end
