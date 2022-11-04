clear
N=225;
c=1;
B=2*pi;
a=1;
n=sqrt(N);


for m=1:n
    H=zeros(N);
B=(m)*2*pi/n;
%פאזה במעבר בציר Y
v=1;
for j=1:n:N-n
for i=1:n
    
    phi_y=((i-1)*a*B/2);
    H(i+j-1,i+j-1+n)=c*exp(-1i*phi_y);
    %תנאי שפה בציר Y
    if j==1 
        H(i+j-1,N-n+i)=c*exp(-1i*(phi_y)/(n-1));
    end
%     if j==1 
%         if rem(v,2)==1
%      H(i+j-1,N-n+i)=c*exp(-1i*phi_y);
%         else  H(i+j-1,N-n+i)=c; 
%         end
%         v=v+1;
%     end
    
end
end 
%פאזה במעבר בציר X
phi_x=0;
w=1;
for j=1:n:N  
for i=1:n-1
    H(i+j-1,i+j)=c*exp(-1i*phi_x);
   
end
%תנאי שפה בציר X
% H(j,j+n-1)=c*exp(-1i*(phi_x)/(n-1));
% if rem(w,2)==1
%  H(j,j+n-1)=c*exp(-1i*phi_x);
% else H(j,j+n-1)=c ;
%end
w=w+1;
  phi_x=phi_x-(B/2*a);
end
H=H+H';
E(:,m)=eig(H);
phi(m)=B*a^2;
end
s=(phi(m)-phi(1))/7;
ph=phi(1):s:phi(m);
figure
plot(E','.k')
xlabel('\phi','fontsize',20)
ylabel('E','fontsize',20)
title('one boundary condition','fontsize',22)
xticklabels(round(ph,2))