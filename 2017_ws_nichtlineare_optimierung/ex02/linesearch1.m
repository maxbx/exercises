function [npoint, nval, ngrad] = linesearch1(fun, cpoint, dir, maxa, cval, cgrad, c1, c2)
%
% performs interpolation search on f
% @fun		... function handle
% cpoint	... curent point
% dir		... descent direction
% maxa		... maximum of alpha
% cval		... current function value
% cgrad		... gradient for the current function value
% c1		... parameter for Armijo condition
% c2		... parameter for curvature condition
%

% it may occur that the stationary point of our interpolation are no Wolfe points for f,
% in that case we can use backtracking, which is commented out, because if there are no
% Wolfe points the program would run forever.

alpha = maxa;
gamma = 0.9;
nval=cval;

%backtracking
while (nval >= (cval+c1*alpha*cgrad'*dir))
  
npoint = cpoint+alpha.*dir;
[nval, ngrad] = (feval(@fun, npoint));

if(ngrad==zeros(length(ngrad),1))
return
end

% calculate coeffitients for cubic interpolation ax^3+bx^2+cx+d
a = 1/(3*alpha^3)*(6*cval+3*alpha*cgrad'*dir-6*nval+3*alpha*ngrad'*dir);
b = 1/alpha^2*(-3*cval-2*alpha*cgrad'*dir+3*nval-alpha*ngrad'*dir);
c = cgrad'*dir;

% calculate stationary points of the interpolation (nur der Minimierer bringt was)

if(a<0) 
alphac = (-b+sqrt(b^2-3*a*c))/(3*a);
elseif(a>0)
alphac = (-b-sqrt(b^2-3*a*c))/(3*a);
elseif((a==0)&&(b<0))
alphac = -c/(2*b);
else
alphac=maxa;  %interpolation yields no useful points for alpha
end

% alpha has to be in [0, maxa]
canidates = [maxa];
if 0 < alphac && alphac < maxa
	canidates = [canidates, alphac];
end

%if 0 < alpha2 && alpha2 < maxa
%	canidates = [canidates, alpha2];
%end

% check the wolfe conditions first for the point with smallest value, then second smallest and biggest
[B, I] = sort(canidates);
for ii = 1:length(I)
	switch I(ii)
	case 1
		alpha = maxa;
	case 2
		alpha = alphac;
	end
	npoint = cpoint+alpha.*dir;
	[nval, ngrad] = (feval(fun, npoint));
	armijo = ( nval <= (cval+c1*alpha*cgrad'*dir) );
	curvature = ( ngrad'*dir >= c2*cgrad'*dir );
	if armijo && curvature
		break;
	end
end
maxa = gamma*maxa;
end
end
