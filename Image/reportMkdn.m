
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%test the QR factorization by testing different A's
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%case 1: decompose an identity matrix
A1=eye(1000,800);
[m1,n1]=size(A1);
[AA1,p1]=house(A1);
[Q1,R1]=getQR(AA1,p1);
nr1=norm(Q1'*Q1-eye(n1));
re1=norm(Q1*R1-A1)/norm(A1);
condA1=cond(A1);


%case 2: decompose a nonsingular matrix
A2=randn(1000,800); %randomly generate a nonsingular matrix
[m2,n2]=size(A2);
[AA2,p2]=house(A2);
[Q2,R2]=getQR(AA2,p2);
nr2=norm(Q2'*Q2-eye(n2));
re2=norm(Q2*R2-A2)/norm(A2);
condA2=cond(A2);


%case 3: decompose an ill-conditioned (singular) matrix 
A3=[1 1 1 1 1; 
    2 2 2 2 2;
    3 3 3 3 3;
    4 4 4 4 4;
    5 5 5 5 5];

[m3,n3]=size(A3);
[AA3,p3]=house(A3);
[Q3,R3]=getQR(AA3,p3);
nr3=norm(Q3'*Q3-eye(n3));
re3=norm(Q3*R3-A3)/norm(A3);
condA3=cond(A3);

%case 4: decompose a random matrix with diagonal entries=0
A4=rand(100,100);
c=diag(A4);
A4=A4-diag(c);
[m4,n4]=size(A4);
[AA4,p4]=house(A4);
[Q4,R4]=getQR(AA4,p4);
nr4=norm(Q4'*Q4-eye(n4));
re4=norm(Q4*R4-A4)/norm(A4);
condA4=cond(A4);

%case 4: decompose a random matrix with upper triangular entries=0
A5=rand(200,100);
B= triu(A5);
A5=A5-B;
[m5,n5]=size(A5);
[AA5,p5]=house(A5);
[Q5,R5]=getQR(AA5,p5);
nr5=norm(Q5'*Q5-eye(n5));
re5=norm(Q5*R5-A5)/norm(A5);
condA5=cond(A5);

%case 5: decompose a random matrix with m>n
A6=rand(1000,800); %randomly generate a nonsingular matrix
[m6,n6]=size(A6);
[AA6,p6]=house(A6);
[Q6,R6]=getQR(AA6,p6);
nr6=norm(Q6'*Q6-eye(n6));
re6=norm(Q6*R6-A6)/norm(A6);
condA6=cond(A6);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%get the computational scaling of the matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%fix k, test n
k=50;
record=zeros(1,20);
for j=1:length(record)
    n=200*j;
    A=rand(n,k);
    f=@() house(A);
    record(j)=timeit(f);
end
x=(1:length(record)).*200;
figure
plot(x,log(record))
title('The computational scaling- change n')
xlabel('n (k=50)') 
ylabel('log(runtime)') 

figure
plot(x,record.^2)
title('The computational scaling- change n')
xlabel('n (k=50)') 
ylabel('time record^2') 

figure
plot(x,record.^(1/5))
title('The computational scaling- change n')
xlabel('n (k=50)') 
ylabel('time record square root') 

%fix n, test k^2
n=600;
record2=zeros(1,20);
for j=1:length(record2)
    k=j*20;
    A=rand(n,k);
    f=@() house(A);
    record2(j)=timeit(f);
end
x2=(1:length(record2)).*20;
figure
%plot(x2,sqrt(record2))
%hold on
%plot(x2,record2.^(1/3))
%hold on
plot(x2,record2)
%hold off
title('The computational scaling- change k^2')
xlabel('k^2 (n=600)') 
ylabel('log(runtime)') 
%legend({'square root','cubic root','original'},'Location','northwest')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Compress the image
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('Cornell_seal_bw.mat')
[B,nr]=comImg(C, 75); %get the compressed image of Cornell seal
imshow(B) %show the compressed image
imshow(C) %show the original image



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Compare the WZ decomposization with SVD decomposition
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
s=svd(C); %get the singular values of matrix C
[U,S,V]=svd(C);
Uk=U(:,1:k);
Vk=V(:,1:k);
Sk=S(1:k,1:k);
D=mat2gray(Uk*Sk*Vk');
imshow(D)

p=min(n1,n2); 
spart=s(k+1:p); %get singular values from k+1 to p
nrsvd= sqrt(sum(spart.^2)); %get the residual norm ||A-SkVkDk||of SVD 
                            %approximation
diff=nr-nrsvd; %get residual norm difference of two methods


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Recover the image
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
betalist=[0.01 0.05 0.1 0.5 0.8 0.99];
%test 1
load('Project1_test1.mat')
for l = 1:length(betalist)
    beta=betalist(l);
    B=recImg(C,75,beta); %get the recovered image
    imwrite(B,strcat('test1_',num2str(beta),'.png')); %save the recovered 
                                                      %image
end

betalist=[0.01 0.05 0.1 0.4 0.7 0.99];
%test 4
load('Project1_test4.mat')
for l = 1:length(betalist)
    beta=betalist(l);
    B=recImg(C,75,beta);
    imwrite(B,strcat('test4_',num2str(beta),'.png'));
end

%test 6
load('Project1_test6.mat')
for l = 1:length(betalist)
    beta=betalist(l);
    B=recImg(C,75,beta);
    imwrite(B,strcat('test6_',num2str(beta),'.png'));
end

