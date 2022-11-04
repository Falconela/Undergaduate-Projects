clear
%   parameters of the problem
DI=0.001;    %  diffusion coefficient of I
KN=2;
KI=0.3;
A=9.5;   % actin mass conservation parameter
DN=0.1;  %  diffusion coefficient of N
xmax=1000;    %   x domain is [0,1000]
tmax=400;    %   t domain is [0,400]
nt=400000;  %   number of steps
dx=xmax/10000;    %   space step
dt=0.01*(dx^2)/(2*DN); %  time step according the stability criterion
nx=xmax;
%   initial matrices
N=zeros(nx+1); 
S=zeros(nx+1);
I=zeros(nx+1);

%   initial conditions
N(1:4,1)=1;
S(:,1)=A;
S(1:4,1)=A-1;
I(1:4,1)=4;

for n=1:nt  %   main loop   
    %   explicit method
    N(2:nx,n+1)=N(2:nx,n)+dt*(DN*((N(3:nx+1,n)+N(1:nx-1,n)-2*N(2:nx,n))/(dx^2))-N(2:nx,n)+(S(2:nx,n).*((N(2:nx,n)).^2))./(1+I(2:nx,n)));     
    S(2:nx,n+1)=S(2:nx,n)+dt*(((S(3:nx+1,n)-2*S(2:nx,n)+S(1:nx-1,n))/(dx^2))+N(2:nx,n)-(S(2:nx,n).*((N(2:nx,n)).^2))./(1+I(2:nx,n))); 
    I(2:nx,n+1)=I(2:nx,n)+dt*(KN*N(2:nx,n)-KI*I(2:nx,n)+DI*((I(3:nx+1,n)-2*I(2:nx,n)+I(1:nx-1,n))/(dx^2)));
    
    %   periodic boundary conditions
    N(nx+1,n+1)=N(2,n+1);
    S(nx+1,n+1)=S(2,n+1);
    I(nx+1,n+1)=I(2,n+1);
    N(1,n+1)=N(nx,n+1);
    S(1,n+1)=S(nx,n+1);
    I(1,n+1)=I(nx,n+1);

%     %   neumann boundary conditions
%     N(1,n+1)=N(2,n+1);
%     N(nx+1,n+1)=N(nx,n+1);
%     S(1,n+1)=S(2,n+1);
%     S(nx+1,n+1)=S(nx,n+1);
%     I(1,n+1)=I(2,n+1);
%     I(nx+1,n+1)=I(nx,n+1);
end
% N(:,1)=N(:,100000);
% S(:,1)=S(:,100000);
% I(:,1)=I(:,100000);
% N(:,2:nt+1)=0;
% S(:,2:nt+1)=0;
% I(:,2:nt+1)=0;
% c=0.004;
% for n=1:nt  %   main loop   
%     %   explicit method
%     N(2:nx,n+1)=N(2:nx,n)+dt*(DN*((N(3:nx+1,n)+N(1:nx-1,n)-2*N(2:nx,n))/(dx^2))-N(2:nx,n)+(S(2:nx,n).*((N(2:nx,n)).^2))./(1+I(2:nx,n)));     
%     S(2:nx,n+1)=S(2:nx,n)+dt*(((S(3:nx+1,n)-2*S(2:nx,n)+S(1:nx-1,n))/(dx^2))+N(2:nx,n)-(S(2:nx,n).*((N(2:nx,n)).^2))./(1+I(2:nx,n))); 
%     I(2:nx,n+1)=I(2:nx,n)+dt*(KN*N(2:nx,n)-KI*I(2:nx,n)+DI*((I(3:nx+1,n)-2*I(2:nx,n)+I(1:nx-1,n))/(dx^2)));
%     
%     %   periodic boundary conditions
%     N(nx+1,n+1)=N(2,n+1);
%     S(nx+1,n+1)=S(2,n+1);
%     I(nx+1,n+1)=I(2,n+1);
%     N(1,n+1)=N(nx,n+1);
%     S(1,n+1)=S(nx,n+1);
%     I(1,n+1)=I(nx,n+1);
% 
%     %   neumann boundary conditions
% %     N(1,n+1)=N(2,n+1);
% %     N(nx+1,n+1)=N(nx,n+1);
% %     S(1,n+1)=S(2,n+1);
% %     S(nx+1,n+1)=S(nx,n+1);
% %     I(1,n+1)=I(2,n+1);
% %     I(nx+1,n+1)=I(nx,n+1);
% end
%  graphical presentation of the amplitude 
% imagesc(N')
% c=gray;
% c=flipud(c);
% colormap(c);
% colorbar
% xlabel('x')
% ylabel('t')