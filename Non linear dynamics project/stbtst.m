%   building the jacobian of the problem
 Jder=zeros(3*nx+3,3*nx+3);
 J0=zeros(3*nx+3,3*nx+3);
 
 Jder(1,1:2)=[(DN/(dxi^2))*-2-(c/dxi) (DN/(dxi^2))+(c/dxi)];
 Jder(2,1:3)=[(DN/(dxi^2)) (DN/(dxi^2))*-2-(c/dxi) (DN/(dxi^2))+(c/dxi)];
 Jder(nx+1,nx:nx+1)=[(DN/(dxi^2)) (DN/(dxi^2))*-2-(c/dxi)];
 J0(1,1)=2*N(1,150000)*S(1,150000)/(1+I(1,150000))-1;
 J0(2,2)=2*N(2,150000)*S(2,150000)/(1+I(2,150000))-1;
 J0(nx+1,nx+1)=2*N(nx+1,150000)*S(nx+1,150000)/(1+I(nx+1,150000))-1;
 for i=3:nx
     Jder(i,i-1:i+1)=Jder(i-1,i-2:i);
     J0(i,i)=2*N(i,150000)*S(i,150000)/(1+I(i,150000))-1;
 end
 Jder(nx+2,nx+2:nx+3)=[(1/(dxi^2))*-2-(c/dxi) (1/(dxi^2))+(c/dxi)];
 Jder(nx+3,nx+2:nx+4)=[(1/(dxi^2)) (1/(dxi^2))*-2-(c/dxi) (1/(dxi^2))+(c/dxi)];
 Jder(2*nx+2,2*nx+1:2*nx+2)=[(1/(dxi^2)) (1/(dxi^2))*-2-(c/dxi)];
 J0(nx+2,nx+2)=2*N(nx+2,150000)*S(nx+2,150000)/(1+I(nx+2,150000))-1;
 J0(nx+3,nx+3)=2*N(nx+3,150000)*S(nx+3,150000)/(1+I(nx+3,150000))-1;
 J0(2*nx+2,2*nx+2)=2*N(2*nx+2,150000)*S(2*nx+2,150000)/(1+I(2*nx+2,150000))-1;
 for i=nx+4:2*nx+1
     Jder(i,i-1:i+1)=Jder(i-1,i-2:i);
     J0(i,i)=2*N(i,150000)*S(i,150000)/(1+I(i,150000))-1;
 end
 Jder(2*nx+3,2*nx+3:2*nx+4)=[(DI/(dxi^2))*-2-(c/dxi) (DI/(dxi^2))+(c/dxi)];
 Jder(2*nx+4,2*nx+3:2*nx+5)=[(DI/(dxi^2)) (DI/(dxi^2))*-2-(c/dxi) (DI/(dxi^2))+(c/dxi)];
 Jder(3*nx+3,3*nx+2:3*nx+3)=[(DI/(dxi^2)) (DI/(dxi^2))*-2-(c/dxi)];
 J0(2*nx+3,2*nx+3)=2*N(2*nx+3,150000)*S(2*nx+3,150000)/(1+I(2*nx+3,150000))-1;
 J0(2*nx+4,2*nx+4)=2*N(2*nx+4,150000)*S(2*nx+4,150000)/(1+I(2*nx+4,150000))-1;
 J0(3*nx+3,3*nx+3)=2*N(3*nx+3,150000)*S(3*nx+3,150000)/(1+I(3*nx+3,150000))-1;
 for i=2*nx+5:3*nx+2
     Jder(i,i-1:i+1)=Jder(i-1,i-2:i);
     J0(i,i)=2*N(i,150000)*S(i,150000)/(1+I(i,150000))-1;
 end
 

    

% J=[2.*N.*S./(1+I)-1 (N.^2)./(1+I) -(S.*N.^2)./((1+I).^2); -2.*N.*S./(1+I)+1 -(N.^2)./(1+I) (S.*N.^2)./((1+I).^2); KN.*ones(nx+1,nx+1) zeros(nx+1,nx+1) -KI.*ones(nx+1,nx+1)];
%   finds the eigenvalues of the jacobian for examination of the stability
% e=eig(J);