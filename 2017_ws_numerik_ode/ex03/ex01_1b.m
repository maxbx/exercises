function ex01_1b

f = @(t, y) t^2/4+y; y0 = 0;
sol = @(t) 1/4*(-2+2*exp(t)-2*t-t.^2);
t0 = 0; T = 1;

N = 12;
maxerror = zeros(N, 3);
times = zeros(N, 3);
for verfahren = 1:3
	switch verfahren
	case 1 % expliziter Euler
		a = 0; c = 1; B = [0];
	case 2 % verbesserter Euler
		a = [0 .5]'; c = [0 1]'; B = [0 0; .5 0];
	case 3 % Verfahren von Heun
		a = [0 1]'; c = [.5 .5]'; B = [0 0; 1 0];
	end

	for n = 0:N
		h = 0.1*2^-n;
		y = erkv(f, t0, T, h, y0, a, c, B);
		tic;
		maxerror(n+1, verfahren) = max(abs(y-sol(t0:h:T)));
		times(n+1, verfahren) = toc;
	end
end
p = [log(0.1*2.^-(0:N)') ones(N+1, 1)]\log(maxerror);
fprintf('geschätzte Konvergenzordnungen:\n expliziter Euler: \t%d\n verbesserter Euler: \t%d\n Verfahren von Heun: \t%d\n', p(1, 1), p(1, 2), p(1, 3));
figure(1);
loglog(0.1*2.^-(0:N), maxerror, '-+');
legend('expliziter Euler', 'verbesserter Euler', 'Verfahren von Heun');
xlabel('h');
ylabel('maxerror');
figure(2);
loglog(times, maxerror, '-+');
legend('expliziter Euler', 'verbesserter Euler', 'Verfahren von Heun');
xlabel('times');
ylabel('maxerror');

end
