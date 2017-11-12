function[xsol,cval,cgrad,it]=BFGS(function_handle,xstart,eps)
%INPUT
%function_handle
%xstart --- starting point
%eps   ---  tolerance for when to stop the iteration

%OUTPUT
%xsol   ---  approximation for a stationary point which fulfills norm(grad(xsol))<eps
%cval,cgrad --- function value / gradient in xsol
%it  ---  number of iterations

[cval,cgrad]=feval(function_handle,xstart);

%choose the identity as a first approximate for the inverse hessian
H_k=eye(length(xstart));
xsol=xstart;
it=0;


while((norm(cgrad))>eps)
%descent direction
p_k=-1*H_k*cgrad;

%compute step length
alpha=linesearchBFGS(function_handle,xsol,p_k,1,cval,cgrad,10^(-4),0.9);
if(alpha==0) %no further progress possible
return
end

xsol=xsol+alpha*p_k;

[nval,ngrad]=feval(function_handle,xsol);


s=alpha*p_k;
y=ngrad-cgrad;

if(it==0)
%Nocedal & Wright's recommendation
H_k=(y'*s)/(y'*y);     
end

%BFGS update
H_k=H_k+(1+(y'*H_k*y)/(s'*y))*((s*s')/(s'*y))-((s*y'*H_k)+(H_k*y*s'))/(s'*y);

cgrad=ngrad;
cval=nval;

%number of iterations
it=it+1;

end

end