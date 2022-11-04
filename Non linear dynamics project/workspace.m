clear
%   parameters of the problem
DI=0.001;    %  diffusion coefficient of I
KN=2;
KI=0.3;
A=9.5;   % actin mass conservation parameter
DN=0.1;  %  diffusion coefficient of N
xmax=1000;    %   x domain is [0,100]
t=4000;    %   t domain is [0,400]
t1=200;%   number of steps
t2=500;
dx=0.1;    %   space step
dt=0.01*(dx^2)/(2*DN); %  time step according the stability criterion
nx=xmax;
%   initial matrices
Nt=zeros(nx+1,t); 
St=zeros(nx+1,t);
It=zeros(nx+1,t);

%  center initial conditions
% Nt(nx/2-2:nx/2+2,1)=1;
% St(:,1)=A;
% St(nx/2-2:nx/2+2,1)=A-1;
% It(nx/2-2:nx/2+2,1)=4;

%  side initial conditions
Nt(1:4,1)=1;
St(:,1)=A;
St(1:4,1)=A-1;
It(1:4,1)=4;
for i=1:t  %   main loop   
    N=zeros(nx+1,t1);
    S=zeros(nx+1,t1);
    I=zeros(nx+1,t1);
    N(:,1)=Nt(:,i);
    S(:,1)=St(:,i);
    I(:,1)=It(:,i);
    %   explicit method
    for n=1:t1
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

    %   neumann boundary conditions
%     N(1,n+1)=N(2,n+1);
%     N(nx+1,n+1)=N(nx,n+1);
%     S(1,n+1)=S(2,n+1);
%     S(nx+1,n+1)=S(nx,n+1);
%     I(1,n+1)=I(2,n+1);
%     I(nx+1,n+1)=I(nx,n+1);
    end
Nt(:,i+1)=N(:,t1+1);
St(:,i+1)=S(:,t1+1);
It(:,i+1)=I(:,t1+1);
clear N S I;
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
% imagesc(Nt')
% c=gray;
% c=flipud(c);
% colormap(c);
% colorbar
% xlabel('x')
% ylabel('t')

h=sum(St)+sum(Nt); %  mass conservation in each step =A
%  finding the velocity of the wave:

% center
% [m1,i1]=max(Nt(nx/2:nx,100));
% [m2,i2]=max(Nt(nx/2:nx,200));
% side
[m1,i1]=max(Nt(:,100));
[m2,i2]=max(Nt(:,200));

c=(-(i2-i1))/(100);
c=0.004*dx/dt;
dxi=dx-c*dt;

%   building the jacobian of the problem
 Jder=zeros(3*nx+3,3*nx+3);
 J0=zeros(3*nx+3,3*nx+3);
 Jder(1,1:2)=[(DN/(dxi^2))*(-2) (DN/(dxi^2))+(c/(2*dxi))];
 Jder(2,1:3)=[(DN/(dxi^2))-(c/(2*dxi)) (DN/(dxi^2))*(-2) (DN/(dxi^2))+(c/(2*dxi))];
 Jder(nx+1,nx:nx+1)=[(DN/(dxi^2))-(c/(2*dxi)) (DN/(dxi^2))*(-2)];
 for i=3:nx
     Jder(i,i-1:i+1)=Jder(i-1,i-2:i);
 end
 Jder(nx+2,nx+2:nx+3)=[(1/(dxi^2))*(-2) (1/(dxi^2))+(c/(2*dxi))];
 Jder(nx+3,nx+2:nx+4)=[(1/(dxi^2))-(c/(2*dxi)) (1/(dxi^2))*(-2) (1/(dxi^2))+(c/(2*dxi))];
 Jder(2*nx+2,2*nx+1:2*nx+2)=[(1/(dxi^2))-(c/(2*dxi)) (1/(dxi^2))*(-2)];
 for i=nx+4:2*nx+1
     Jder(i,i-1:i+1)=Jder(i-1,i-2:i);
 end
 Jder(2*nx+3,2*nx+3:2*nx+4)=[(DI/(dxi^2))*(-2) (DI/(dxi^2))+(c/(2*dxi))];
 Jder(2*nx+4,2*nx+3:2*nx+5)=[(DI/(dxi^2))-(c/(2*dxi)) (DI/(dxi^2))*(-2) (DI/(dxi^2))+(c/(2*dxi))];
 Jder(3*nx+3,3*nx+2:3*nx+3)=[(DI/(dxi^2))-(c/(2*dxi)) (DI/(dxi^2))*(-2)];
 for i=2*nx+5:3*nx+2
     Jder(i,i-1:i+1)=Jder(i-1,i-2:i);
 end
 for i=1:nx+1
     J0(i,i)=(2*Nt(i,t2)*St(i,t2))/(1+It(i,t2))-1;
     J0(nx+1+i,nx+1+i)=-((Nt(i,t2))^2)/(1+It(i,t2));
     J0(2*nx+2+i,2*nx+2+i)=-KI;
     J0(i,nx+1+i)=((Nt(i,t2))^2)/(1+It(i,t2));
     J0(i,2*nx+2+i)=-(St(i,t2)*((Nt(i,t2))^2)/((1+It(i,t2))^2));
     J0(nx+1+i,i)=-(2*Nt(i,t2)*St(i,t2))/(1+It(i,t2))+1;
     J0(2*nx+2+i,i)=KN;
     J0(nx+1+i,2*nx+2+i)=(St(i,t2)*(Nt(i,t2)^2)/((1+It(i,t2))^2));
 end

 
    %  neumann boundary conditions for the jacobian
%  J(1,1:3*nx+3)=J(2,1:3*nx+3);
%  J(nx+1,1:3*nx+3)=J(nx,1:3*nx+3);
%  J(nx+2,1:3*nx+3)=J(nx+3,1:3*nx+3);
%  J(2*nx+2,1:3*nx+3)=J(2*nx+1,1:3*nx+3);
%  J(2*nx+3,1:3*nx+3)=J(2*nx+4,1:3*nx+3);
%  J(3*nx+3,1:3*nx+3)=J(3*nx+2,1:3*nx+3);

    %  periodic boundary conditions for the jacobian
 J0(nx+1,1:3*nx+3)=J0(2,1:3*nx+3);
 J0(2*nx+2,1:3*nx+3)=J0(nx+3,1:3*nx+3);
 J0(3*nx+3,1:3*nx+3)=J0(2*nx+4,1:3*nx+3);
 J0(1,1:3*nx+3)=J0(nx,1:3*nx+3);
 J0(nx+2,1:3*nx+3)=J0(2*nx+1,1:3*nx+3);
 J0(2*nx+3,1:3*nx+3)=J0(3*nx+2,1:3*nx+3);
 Jder(1,:)=0;Jder(2,:)=0;Jder(nx,:)=0;Jder(nx+1,:)=0;Jder(nx+2,:)=0;Jder(nx+3,:)=0;
 Jder(2*nx+1,:)=0;Jder(2*nx+2,:)=0;Jder(2*nx+3,:)=0;Jder(2*nx+4,:)=0;Jder(3*nx+2,:)=0;Jder(3*nx+3,:)=0;
 Jder(1,2)=DN/(dxi^2)+c/(2*dxi);Jder(1,nx)=-2*DN/(dxi^2);
 Jder(nx+2,nx+3)=1/(dxi^2)+c/(2*dxi);Jder(nx+2,2*nx+1)=-2/(dxi^2);
 Jder(2*nx+3,2*nx+4)=DI/(dxi^2)+c/(2*dxi);Jder(2*nx+3,3*nx+2)=-2*DI/(dxi^2);
 Jder(2,2)=-2*DN/(dxi^2);Jder(2,3)=DN/(dxi^2)+c/(2*dxi);Jder(2,nx)=DN/(dxi^2)-c/(2*dxi);
 Jder(nx+3,nx+3)=-2/(dxi^2);Jder(nx+3,nx+4)=1/(dxi^2)+c/(2*dxi);Jder(nx+3,2*nx+1)=1/(dxi^2)-c/(2*dxi);
 Jder(2*nx+4,2*nx+4)=-2*DI/(dxi^2);Jder(2*nx+4,2*nx+5)=DI/(dxi^2)+c/(2*dxi);Jder(2*nx+4,3*nx+2)=DI/(dxi^2)-c/(2*dxi);
 Jder(nx,2)=DN/(dxi^2)+c/(2*dxi);Jder(nx,nx-1)=DN/(dxi^2)-c/(2*dxi);Jder(nx,nx)=-2*DN/(dxi^2);
 Jder(2*nx+1,nx+3)=1/(dxi^2)+c/(2*dxi);Jder(2*nx+1,2*nx)=1/(dxi^2)-c/(2*dxi);Jder(2*nx+1,2*nx+1)=-2*1/(dxi^2);
 Jder(3*nx+2,2*nx+4)=DI/(dxi^2)+c/(2*dxi);Jder(3*nx+2,3*nx+1)=DI/(dxi^2)-c/(2*dxi);Jder(3*nx+2,3*nx+2)=-2*DI/(dxi^2);
 Jder(nx+1,2)=-2*DN/(dxi^2);Jder(nx+1,nx)=DN/(dxi^2)-c/(2*dxi);
 Jder(2*nx+2,nx+3)=-2*1/(dxi^2);Jder(2*nx+2,2*nx+1)=1/(dxi^2)-c/(2*dxi);
 Jder(3*nx+3,2*nx+4)=-2*DI/(dxi^2);Jder(3*nx+3,3*nx+2)=DI/(dxi^2)-c/(2*dxi);
 J=Jder+J0;
[d,e]=eigs(J,5,'lr');
x=sort((eig(J)));
plot(1:1001,d(1:1001,1),1:1001,Nt(:,t2))
