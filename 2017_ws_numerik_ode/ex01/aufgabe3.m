function aufgabe3

alpha1 = 1; alpha2 = 2; beta1 = 0.05; beta2 = 0.05;
f = @(x, y) [alpha1*y(1)-beta1*y(1)*y(2), -alpha2*y(2)+beta2*y(1)*y(2)];
y0 = [100 10];
t0 = 0; T = 20; h = 0.001;

y = ExpliziterEuler(f, t0, T, h, y0);
plot(linspace(t0, T, size(y, 1)), y(:, :)', 'k');
set(gca, 'FontSize', 20);
xlabel('Zeit t', 'FontSize', 20);
ylabel('Population y', 'FontSize', 20);

end
