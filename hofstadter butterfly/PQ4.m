clear
n=9; %number of sites
c=1; %transition amplitude
a=1; %distance between two adjacent sites 
N=n^2;

for m=1:n
    H=zeros(N);
 B=((m)*2*pi/((n))); % magnetic field
%transition in y axis
for j=1:n:N-n
for i=1:n
    phi_y=((i-1)*a*B);
    H(i+j-1,i+j-1+n)=c*exp(1i*phi_y);
%     PBC in y axis
    if j==1 
        H(i+j-1,N-n+i)=c*exp(-1i*(phi_y));
    end
end
end 
%transition in x axis
for j=1:n:N  
for i=1:n-1
    H(i+j-1,i+j)=c;  
end
%PBC in x axis
H(j,j+n-1)=c;
end
H=H+H';
E(:,m)=eig(H);
phi(m)=B*a^2;
end

figure
cc=hsv(size(E,1));
for q=1:size(E,1)
    hold on
 plot(E(q,:),'.','color',cc(q,:),'LineWidth',2)
end
 
s=(phi(m)-phi(1))/8;
ph=phi(1):s:phi(m);
xlabel('\Phi','fontsize',20)
ylabel('Energy','fontsize',20)
title('PBC in the x axis & OBC in the y axis','fontsize',22)
xticklabels(round(ph,2))