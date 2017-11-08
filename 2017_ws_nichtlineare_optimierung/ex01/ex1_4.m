function rosenbrock(x1,y1);
x=-1:0.02:2;        % grid points along x-axes         
y=-0.5:0.02:3;      % grid points along y-axes

nx=length(x);       % size of grid along x-axes 
ny=length(y);       % siye of grid along y-axes

%
% evaluate Rosenbrock function in grid points
%

g=[2*x1-2-400*(y1-x1^2)*x1;200*(y1-x1^2)];
h=[2+1200*x1^2-400*y1, -400*x1;-400*x1,200];
f=[100*(y1-x1^2)^2+(1-x1)^2];


fz=zeros(nx,ny);    % intiliaze matrix size to grid size and set to zero  
lm=zeros(nx,ny);    % intiliaze matrix size to grid size and set to zero  
qm=zeros(nx,ny);    % intiliaze matrix size to grid size and set to zero  
U=zeros(nx,ny);     % intiliaze matrix size to grid size and set to zero  
V=zeros(nx,ny);     % intiliaze matrix size to grid size and set to zero  
for i=1:nx
  for j=1:ny
    % compute f(x(i),y(j)) and store in fz(i,j)
    fz(i,j)= 100*(y(j)-x(i)^2)^2+(1-x(i))^2;
    U(i,j)=2*x(i)-2-400*(y(j)-x(i)^2)*x(i);
    V(i,j)=200*(y(j)-x(i)^2);
    lm(i,j)=f+g'*[x(i)-x1;y(j)-y1];
    qm(i,j)=f+g'*[x(i)-x1;y(j)-y1]+[x(i)-x1;y(j)-y1]'*h*[x(i)-x1;y(j)-y1]/2;
  end
end

%
% make a contour plot for prespecified contour levels
%

figure(2);               % output next plot on window of figure 1 
hold off;                % overwrite old plot in figure 1
contour(x,y,fz',[0.1,0.5,2.5,12.5,62.5,312.5]);     % plot the function
xlabel('x');             % output 'x' along x axis
ylabel('y');             % output 'y' along y axis
hold on;
plot([x1],[y1],'or');

figure(1);
clf;
hold off;
axis([min(x) max(x) min(y) max(y) -200 1000 -200 1000]);
view(3);
hold on;

I=1:ceil(length(x)/30):length(x);
J=1:ceil(length(y)/30):length(y);
mesh(x(I),y(J),fz(I,J)');

figure(3);
clf;
hold off;
quiver(x(I)'*ones(1,length(J)),ones(length(I),1)*y(J),U(I,J),V(I,J),5);


pause 


figure(2);
hold on;


contour(x,y,qm',[0.1,0.5,2.5,12.5,62.5,312.5]);     % plot the quadratic model

%
% make a mesh plot
%

figure(1)
hold on;

mesh(x(I),y(J),qm(I,J)');

