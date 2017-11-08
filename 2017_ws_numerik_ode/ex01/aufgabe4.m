function aufgabe4

N = 30;
f = @(x, y) (toeplitz([-1 zeros(1, 2*N-3) 1 0], [-1 0 1 zeros(1, 2*N-3)])*y')';

y0 = rand(N, 2);

t0 = 0; T = 10; h = 0.1;

y = ExpliziterEuler(f, t0, T, h, y0);
for ii = 1:size(y, 1)
	scatter(y(ii, 1:2:end), y(ii, 2:2:end), 100, 'filled');
	axis([0 1, 0 1]);
	axis off;
	drawnow;
%	print('-dpng', '-r50', sprintf('%.5d.png', ii));
end

end
