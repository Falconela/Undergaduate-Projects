%   parameters of the problem
DI=0.001;    %  diffusion coefficient of I
KN=2;
KI=0.3;
A=7.7;   % actin mass conservation parameter
DN=0.1;  %  diffusion coefficient of N
xmax=100    %   x domain is [0,100]
tmax=400    %   t domain is [0,400]
nt=100;  %   number of time steps
dt=tmax/nt;    %   time step
nx=100;  %   number of space steps
dx=xmax/nx;    %   space step
%   initial matrices
N=zeros(nx+1,nt+1); 
S=zeros(nx+1,nt+1);
I=zeros(nx+1,nt+1);

%   initial conditions 0
for  j=1:nx+1
   x(j)=(j-1)*dx;
   N(j,1)=10^-3; 
   S(1,1)=A-10^-3;
   I(j,1)=10^-3;
end

%   initial conditions +
% B=0.5*(A-KN/KI+sqrt((A-KN/KI)^2-4))
% for  j=1:nx+1
%    N(j,1)=B; 
%    S(1,1)=A-B;
%    I(j,1)=KN/KI*B;
% end

%   boundary conditions
for n=1:nt+1
    N(max(j),n)=N(1,n);
    S(max(j),n)=S(1,n);
    I(max(j),n)=I(1,n);
    time(n)=(n-1)*dt;
end

%   explicit method
if 2*DI*dt/(dx*dx)<=1  %   stability=<1
    for n=1:nt  %   time loop
        for j=2:nx  %   space loop- 
            %if sum(sum(N))+sum(sum(S))==A %   int [N(x)+S(x)]=A
                N(j,n+1)=N(j,n)+dt*(DN*((N(j+1,n)-2*N(j,n)+N(j-1,n))/(dx^2))-N(j,n)+(S(n,j)*(N(j,n))^2)/(1+I(j,n)));     
                S(j,n+1)=S(j,n)+dt*(((S(j+1,n)-2*S(j,n)+S(j-1,n))/(dx^2))+N(j,n)-(S(n,j)*(N(j,n))^2)/(1+I(j,n))); 
                I(j,n+1)=I(j,n)+dt*(KN*N(j,n)-KI*I(j,n)+DI*((I(j+1,n)-2*I(j,n)+I(j-1,n))/(dx^2)));
            %end
        end
    end
end


