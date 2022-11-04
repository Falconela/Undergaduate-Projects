clear
N=9;
c=1;
B=1;
a=1;
n=sqrt(N);
%כיול לנדאו ותנאי שפה בציר X
%pi=B(0,x,0)

for m=1:n
    H=zeros(N);
B=pi/4;
% (m*2*pi/((n));

%פאזה במעבר בציר Y
for j=1:n:N-n
for i=1:n
    phi_y=((i-1)*a*B);
    H(i+j-1,i+j-1+n)=c*exp(1i*phi_y);
%     boundary in y axis
    if j==1 
        H(i+j-1,N-n+i)=c*exp(-1i*(phi_y));
    end
end
end 
%מעברים בציר X
for j=1:n:N  
for i=1:n-1
    H(i+j-1,i+j)=c;  
end
%תנאי שפה בציר X
% H(j,j+n-1)=c;
end
H=H+H';
E(:,m)=eig(H);
phi(m)=B*a^2;
end
figure
cc=hsv(length(E));
for q=1:length(E)
    hold on
 plot(E(q,:),'.','color',cc(q,:),'LineWidth',2)
end
% 
s=(phi(m)-phi(1))/7;
ph=phi(1):s:phi(m);
% figure
% plot(E','.k')
xlabel('\phi','fontsize',20)
ylabel('E','fontsize',20)
title('one boundary condition','fontsize',22)
xticklabels(round(ph,2))