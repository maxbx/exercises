function y = erkv(f, t0, T, h, y0, a, c, B)

S = length(a);
N = floor((T-t0)/h)+1;
y = [y0, zeros(length(y0), N-1)];
t = t0;

for ii = 1:N-1
	K = zeros(length(y0), S);
	for jj = 1:S
		K(jj) = f(t+a(jj)*h, y(ii)+h*sum(K*B(jj, :)', 2));
	end
	y(:, ii+1) = y(:, ii)+h*(K*c);
	t = t+h;
end

end
