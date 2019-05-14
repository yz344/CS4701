function [B,nr]=comImg(A,k)
%this function returns the compressed image matrix of A given compression 
%rank k and also returns the residual norm of ||A-WZ'||.
[n1,n2]=size(A);
W=rand(n1,k); %initialize W
ZT=zeros(k,n2); %initialize Z
iter=1;
nr=1;
while iter<50 && nr > 1e-6
    %update Z
    [Q1,p1]=house(W);
    for i=1:n2
         ZT(:,i)=lsol2(A(:,i),Q1,p1);
    end
    
    %update W
    Z=ZT';
    WT=W';
    AT=A';
    [Q2,p2]=house(Z);
    for j=1:n1
         [WT(:,j),nr]=lsol(AT(:,j),Q2,p2); 
    end
    W=WT';
    iter=iter+1;

end

nr=norm(A-W*Z','fro');
B=mat2gray(W*Z');
end

    
