function y = ExpliziterEuler(f, t0, T, h, y0)
	t = t0;

	ii = 1;
	y(ii, :) = reshape(y0, 1, []);
	
	while t <= T
		y(ii+1, :) = y(ii, :)+h*f(t, y(ii, :));
		ii = ii+1;
		t = t+h;
	end

end
