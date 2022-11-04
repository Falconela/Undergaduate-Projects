 clear
%parameters
N=9;
c=1;
a=1;
for m=1:100
B=(0.01*m-0.01);
% constract H for no bwondery condithion
H=diag(ones(N-1,1),1)*c;
H=H+diag(ones(N-sqrt(N),1),sqrt(N))*c;
H=H+diag(ones(sqrt(N),1),-N+sqrt(N))*c;
for j=1:N
    for k=0:sqrt(N)-1
        for l=1:sqrt(N)
            phi_x=-B/2*k;
            H(k*sqrt(N)+l,j)= H(k*sqrt(N)+l,j)*exp(-1i*phi_x);
        end
    end
end
% for l=1:N
%     for k=0:sqrt(N)-1
%         for j=1:sqrt(N)
%             phi_y=B/2*k;
%             H(l,k*sqrt(N)+j)= H(l,k*sqrt(N)+j)*exp(-1i*phi_y);
%         end
%     end
% end
H=H+H';
E(:,m)=eig(H);
end
B=0:224;
plot(E,B,'.')



%constract H with 1 bwondery condition
% for m=1:100
%  B=0.01*m-0.01;
% h=diag(ones(N-1,1),1)*c;
% %%h=h+diag(ones(N-1,1),-1)*c;
% %%h=h+diag(ones(N-sqrt(N),1),-sqrt(N))*c;
% h=h+diag(ones(N-sqrt(N),1),sqrt(N))*c;
% h=h+diag(ones(sqrt(N),1),-N+sqrt(N))*c;
% %h=h+diag(ones(sqrt(N),1),N-sqrt(N))*c;
% for j=1:N
%     for l=1:N
%         phi_x=-B/2*j*a;
%         h(l,j)= h(l,j)*exp(-i*phi_x);
%     end
% 
% end
% for l=1:N
%     for j=1:N
%         phi_y=B/2*l*a;
%         h(l,j)=h(l,j)*exp(-i*phi_y);
%     end
% end
% E(:,m)=eig(h);
% h=h+h';
% end
% plot(E,'.k')
