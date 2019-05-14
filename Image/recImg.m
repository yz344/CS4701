function B=recImg(A,mask,k,beta)
%this function returns a recovered image given an incomplete image matrix
%A and the mask matrix, rank k and l2 norm parameter beta.
[n1,n2]=size(A);
nullmat=eye(k)*beta; %prepare several matrices for least square problem
nullmat1=ones(k,n2);
nullmat2=ones(k,n1);
nullvec=zeros(k,1);

W=rand(n1,k); %initialize W
M=[W;nullmat]; %initialize M used to solve Mz=b
ZT=zeros(k,n2); %initialize Z
iter=1;
newMask=[mask;nullmat1];
newMaskT=[mask';nullmat2];
while iter<=100 
    %update Z   
    for i=1:n2
        b=[A(:,i);nullvec];
        rl=find(newMask(:,i)==1); %get a rowlist where A is not equal to 1
        [Q1,p1]=house(M(rl,:));
        ZT(:,i)=lsol2(b(rl),Q1,p1); 
        %get the solution from selected pixels
        
    end
    
    %update W
    
    Z=ZT';
    WT=W';
    N=[Z;nullmat]; %initialize N used to solve Nz=b2
    AT=A';
    
    for j=1:n1
         b2=[AT(:,j);nullvec];
         rl2=find(newMaskT(:,j)==1); %get a rowlist where A is not equal to 1
         [Q2,p2]=house(N(rl2,:));
         WT(:,j)=lsol2(b2(rl2),Q2,p2); 
    end
    W=WT';
    M=[W;nullmat];
    iter=iter+1;

end

Z=N(1:n2,:);
B=mat2gray(W*Z');

end

