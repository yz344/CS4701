function [Q,R]=getQR(AA,p) 
%get Q and R from house function

[m,n]=size(AA); %suppose m>n
R=triu(AA(1:n,:));
Q=eye(m);
for k=1:n
    w=zeros(m,1);
    
    w(k:m,1) = [p(k);AA(k+1:m,k)];
    Q=Q-2*Q*w*w'; 

end

Q=Q(:,1:n);
end