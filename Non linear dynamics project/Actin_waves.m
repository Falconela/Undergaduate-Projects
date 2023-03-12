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
N=zeros(nx+1,nt+1); 
S=zeros(nx+1,nt+1);
I=zeros(nx+1,nt+1);

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
% xlabel('x')
% ylabel('t')

h=sum(S)+sum(N); %  mass conservation in each step =A
%  finding the velocity of the wave
[m1,i1]=max(N(:,50000));
[m2,i2]=max(N(:,100000));
c=((i2-i1)*dx)/(50000*dt);
% for ni=1:20
% c=ni*0.05+0.45;
% c=0.9;
dxi=dx-c*dt;

%   building the jacobian of the problem
 Jder=zeros(3*nx+3,3*nx+3);
 J0=zeros(3*nx+3,3*nx+3);
 Jder(1,1:2)=[(DN/(dxi^2))*(-2) (DN/(dxi^2))+(c/(2*dxi))];
  Jder(1,nx+1)=(DN/(dxi^2))-(c/(2*dxi));
 Jder(2,1:3)=[(DN/(dxi^2))-(c/(2*dxi)) (DN/(dxi^2))*(-2) (DN/(dxi^2))+(c/(2*dxi))];
  Jder(nx+1,1)=(DN/(dxi^2))+(c/(2*dxi));
 Jder(nx+1,nx:nx+1)=[(DN/(dxi^2))-(c/(2*dxi)) (DN/(dxi^2))*(-2)];

 for i=3:nx
     Jder(i,i-1:i+1)=Jder(i-1,i-2:i);
 end
 Jder(nx+2,nx+2:nx+3)=[(1/(dxi^2))*(-2) (1/(dxi^2))+(c/(2*dxi))];
  Jder(nx+2,2*nx+2)=(1/(dxi^2))-(c/(2*dxi));
 Jder(nx+3,nx+2:nx+4)=[(1/(dxi^2))-(c/(2*dxi)) (1/(dxi^2))*(-2) (1/(dxi^2))+(c/(2*dxi))];
  Jder(2*nx+2,nx+2)=(1/(dxi^2))+(c/(2*dxi));
 Jder(2*nx+2,2*nx+1:2*nx+2)=[(1/(dxi^2))-(c/(2*dxi)) (1/(dxi^2))*(-2)];
 for i=nx+4:2*nx+1
     Jder(i,i-1:i+1)=Jder(i-1,i-2:i);
 end
 Jder(2*nx+3,2*nx+3:2*nx+4)=[(DI/(dxi^2))*(-2) (DI/(dxi^2))+(c/(2*dxi))];
  Jder(2*nx+3,3*nx+3)=(DI/(dxi^2))-(c/(2*dxi));
 Jder(2*nx+4,2*nx+3:2*nx+5)=[(DI/(dxi^2))-(c/(2*dxi)) (DI/(dxi^2))*(-2) (DI/(dxi^2))+(c/(2*dxi))];
  Jder(3*nx+3,2*nx+3)=(DI/(dxi^2))+(c/(2*dxi));
 Jder(3*nx+3,3*nx+2:3*nx+3)=[(DI/(dxi^2))-(c/(2*dxi)) (DI/(dxi^2))*(-2)];
 for i=2*nx+5:3*nx+2
     Jder(i,i-1:i+1)=Jder(i-1,i-2:i);
 end
 t=50000;
 N(500:1001,t)=N(500,t);
 S(500:1001,t)=S(500,t);
 I(500:1001,t)=I(500,t);
 for i=1:nx+1
     J0(i,i)=(2*N(i,t)*S(i,t))/(1+I(i,t))-1;
     J0(nx+1+i,nx+1+i)=-((N(i,t))^2)/(1+I(i,t));
     J0(2*nx+2+i,2*nx+2+i)=-KI;
     J0(i,nx+1+i)=((N(i,t))^2)/(1+I(i,t));
     J0(i,2*nx+2+i)=-(S(i,t)*((N(i,t))^2)/((1+I(i,t))^2));
     J0(nx+1+i,i)=-(2*N(i,t)*S(i,t))/(1+I(i,t))+1;
     J0(2*nx+2+i,i)=KN;
     J0(nx+1+i,2*nx+2+i)=(S(i,t)*(N(i,t)^2)/((1+I(i,t))^2));
 end

 
    %  neumann boundary conditions for the jacobian
%  J(1,1:3*nx+3)=J(2,1:3*nx+3);
%  J(nx+1,1:3*nx+3)=J(nx,1:3*nx+3);
%  J(nx+2,1:3*nx+3)=J(nx+3,1:3*nx+3);
%  J(2*nx+2,1:3*nx+3)=J(2*nx+1,1:3*nx+3);
%  J(2*nx+3,1:3*nx+3)=J(2*nx+4,1:3*nx+3);
%  J(3*nx+3,1:3*nx+3)=J(3*nx+2,1:3*nx+3);

    %  periodic boundary conditions for the jacobian
%  J0(nx+1,:)=J0(2,:);
%  J0(2*nx+2,:)=J0(nx+3,:);
%  J0(3*nx+3,:)=J0(2*nx+4,:);
%  J0(1,:)=J0(nx,:);
%  J0(nx+2,:)=J0(2*nx+1,:);
%  J0(2*nx+3,:)=J0(3*nx+2,:);
%  Jder(nx+1,:)=Jder(2,:);
%  Jder(2*nx+2,:)=Jder(nx+3,:);
%  Jder(3*nx+3,:)=Jder(2*nx+4,:);
%  Jder(1,:)=Jder(nx,:);
%  Jder(nx+2,:)=Jder(2*nx+1,:);
%  Jder(2*nx+3,:)=Jder(3*nx+2,:);
%  Jder(1,:)=0;Jder(2,:)=0;Jder(nx,:)=0;Jder(nx+1,:)=0;Jder(nx+2,:)=0;Jder(nx+3,:)=0;
%  Jder(2*nx+1,:)=0;Jder(2*nx+2,:)=0;Jder(2*nx+3,:)=0;Jder(2*nx+4,:)=0;Jder(3*nx+2,:)=0;Jder(3*nx+3,:)=0;
%  Jder(1,2)=DN/(dxi^2)+c/(2*dxi);Jder(1,nx)=-2*DN/(dxi^2);
%  Jder(nx+2,nx+3)=1/(dxi^2)+c/(2*dxi);Jder(nx+2,2*nx+1)=-2/(dxi^2);
%  Jder(2*nx+3,2*nx+4)=DI/(dxi^2)+c/(2*dxi);Jder(2*nx+3,3*nx+2)=-2*DI/(dxi^2);
%  Jder(2,2)=-2*DN/(dxi^2);Jder(2,3)=DN/(dxi^2)+c/(2*dxi);Jder(2,nx)=DN/(dxi^2)-c/(2*dxi);
%  Jder(nx+3,nx+3)=-2/(dxi^2);Jder(nx+3,nx+4)=1/(dxi^2)+c/(2*dxi);Jder(nx+3,2*nx+1)=1/(dxi^2)-c/(2*dxi);
%  Jder(2*nx+4,2*nx+4)=-2*DI/(dxi^2);Jder(2*nx+4,2*nx+5)=DI/(dxi^2)+c/(2*dxi);Jder(2*nx+4,3*nx+2)=DI/(dxi^2)-c/(2*dxi);
%  Jder(nx,2)=DN/(dxi^2)+c/(2*dxi);Jder(nx,nx-1)=DN/(dxi^2)-c/(2*dxi);Jder(nx,nx)=-2*DN/(dxi^2);
%  Jder(2*nx+1,nx+3)=1/(dxi^2)+c/(2*dxi);Jder(2*nx+1,2*nx)=1/(dxi^2)-c/(2*dxi);Jder(2*nx+1,2*nx+1)=-2*1/(dxi^2);
%  Jder(3*nx+2,2*nx+4)=DI/(dxi^2)+c/(2*dxi);Jder(3*nx+2,3*nx+1)=DI/(dxi^2)-c/(2*dxi);Jder(3*nx+2,3*nx+2)=-2*DI/(dxi^2);
%  Jder(nx+1,2)=-2*DN/(dxi^2);Jder(nx+1,nx)=DN/(dxi^2)-c/(2*dxi);
%  Jder(2*nx+2,nx+3)=-2*1/(dxi^2);Jder(2*nx+2,2*nx+1)=1/(dxi^2)-c/(2*dxi);
%  Jder(3*nx+3,2*nx+4)=-2*DI/(dxi^2);Jder(3*nx+3,3*nx+2)=DI/(dxi^2)-c/(2*dxi);
 J=Jder+J0;
%  opts.p=100;
%  opts.maxit=500;
% [d,e]=eigs(J,5,'lr',opts);
 [V,D]=eig(full(J));
 q=real(diag(D));
 [z,ind]=sort(q(1:1001),'descend');
% plot(1:1001,d(1:1001,1),1:1001,N(:,t))
% plot(1:1001,d(1:1001,1),1:1001,d(1002:2002,1),1:1001,d(2003:3003,1))
% legend('N','S','I')
% l(:,ni)=q(1:2);
% end
% Nt2=zeros(1002,1);
% Nt1=zeros(1002,1);
% Nt2(2:1002)=N(:,t);
% Nt1(1:1001)=N(:,t);
% 
% plot(1:1002,(Nt2-Nt1)/dx)
% plot(1:1001,N(:,t),'b',1:1001,S(:,t),'r',1:1001,I(:,t),'g')
plot(1:1001,sort(q(1:1001)),'ko')
