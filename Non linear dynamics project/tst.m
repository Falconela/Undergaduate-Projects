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
     J0(i,i)=(2*N(i,300000)*S(i,300000))/(1+I(i,300000))-1;
     J0(nx+1+i,nx+1+i)=-((N(i,300000))^2)/(1+I(i,300000));
     J0(2*nx+2+i,2*nx+2+i)=-KI;
     J0(i,nx+1+i)=((N(i,300000))^2)/(1+I(i,300000));
     J0(i,2*nx+2+i)=-(S(i,300000)*((N(i,300000))^2)/((1+I(i,300000))^2));
     J0(nx+1+i,i)=-(2*N(i,300000)*S(i,300000))/(1+I(i,300000))+1;
     J0(2*nx+2+i,i)=KN;
     J0(nx+1+i,2*nx+2+i)=(S(i,300000)*(N(i,300000)^2)/((1+I(i,300000))^2));
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