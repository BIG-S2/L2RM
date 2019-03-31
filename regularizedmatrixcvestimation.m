%%% n is the sample size, p, q are the dimensions of the 2D image, s is the
%%% dimension of the selected covariates after screening, not the original
%%% dimension
%%% totalX is n by s, totalY is a np*q matrix.
%%% The (i*p-p+1):(i*p) rows of totalY correspond to the p by q response from the ith subject
function[regularizedB]=regularizedmatrixcvestimation(totalX,totalY,cvsplit,p,q,s,n) 
%tic;

lambdavalue=linspace(-5,5,31);
lengthlambda=length(lambdavalue);
candidateB=[];
testerror=[];

jj=1;
element=[];
element(1)=1;

while jj<=lengthlambda & element(jj)>0
lambda=exp(lambdavalue(jj));

for rr=1:5,
    trainindex=find(cvsplit.training(rr)==1);
    testindex=find(cvsplit.test(rr)==1);
    X=totalX(trainindex,:);
    testX=totalX(testindex,:);
    
    Y=[];
    testY=[];

    for j=1:length(trainindex),
        Y((p*j-p+1):(p*j),:)=totalY((p*trainindex(j)-p+1):(p*trainindex(j)),:);
    end
    
    for j=1:length(testindex),
        testY((p*j-p+1):(p*j),:)=totalY((p*testindex(j)-p+1):(p*testindex(j)),:);
    end

%profile on
alpha=[];
B=[];
B{1}=zeros(p,s*q);
B0=zeros(p,s*q);

tempn=length(trainindex);

%profile on
nrm = (norm(X,2))^2;

delta=tempn/nrm;
alpha0=0;
alpha{1}=1;
solution=[];
temphvalue=[];



t=1;

while t==1 || t==2 || abs(temphvalue{t-1}-temphvalue{t-2})/abs(temphvalue{t-2})>10^(-4),
    if t==1,
       solution{t}=B{t}+(alpha0-1)/alpha{t}*(B{t}-B0);
    end
    if t>1,
       solution{t}=B{t}+(alpha{t-1}-1)/alpha{t}*(B{t}-B{t-1});
    end   
           Btemp=zeros(p,s*q);
           Atemp=solution{t}-delta*dlossmatrixdividen(Y,X,solution{t},s,length(trainindex));
           subAtemp=[];
           subBtemp=[];
           for j=1:s,
               subAtemp(:,:,j)=Atemp(:,(j*q-q+1):(j*q));
               [U,S,V]=svd(subAtemp(:,:,j),0);
               avector=(diag(S));
               bvector=(avector-lambda*delta*ones(length(avector),1)).*(avector>lambda*delta*ones(length(avector),1));
               if p==q,
                  subBtemp(:,:,j)=U*diag(bvector)*V';
               else
                  BS=zeros(p,q);
                  BS(1:min(p,q),1:min(p,q))=diag(bvector);
                  subBtemp(:,:,j)=U*BS*V';
               end
           end
     Btemp=reshape(subBtemp,p,s*q);          
           %hfunction(Y,X,Btemp,lambda)
     temphvalue{t}=hmatrixfunctiondividen(Y,X,B{t},lambda,s,tempn);
     B{t+1}=Btemp;
     alpha{t+1}=(1+sqrt(1+(2*alpha{t})^2))/2; 
     t=t+1;
end  

%estimateB1=B{t}(:,1:q);
%estimateB2=B{t}(:,(q+1):(2*q));

%profile viewer

tempcandidateB=B{t};
temptesterror(rr)=lossmatrixdividen(testY,testX,tempcandidateB,s,tempn);
tempelement(rr)=max(max(abs(tempcandidateB)));
end

testerror(jj)=mean(temptesterror);

jj=jj+1;
element(jj)=max(abs(tempelement));
end

testerrorindex=find(testerror==min(testerror));

lambda=exp(lambdavalue(testerrorindex));

%profile on
alpha=[];
B=[];
B{1}=zeros(p,s*q);
B0=zeros(p,s*q);

%profile on
nrm = (norm(totalX,2))^2;

delta=n/nrm;
alpha0=0;
alpha{1}=1;
solution=[];
temphvalue=[];



t=1;

while t==1 || t==2 || abs(temphvalue{t-1}-temphvalue{t-2})/abs(temphvalue{t-2})>10^(-4),
    if t==1,
       solution{t}=B{t}+(alpha0-1)/alpha{t}*(B{t}-B0);
    end
    if t>1,
       solution{t}=B{t}+(alpha{t-1}-1)/alpha{t}*(B{t}-B{t-1});
    end   
           Btemp=zeros(p,s*q);
           Atemp=solution{t}-delta*dlossmatrixdividen(totalY,totalX,solution{t},s,n);
           subAtemp=[];
           subBtemp=[];
           for j=1:s,
               subAtemp(:,:,j)=Atemp(:,(j*q-q+1):(j*q));
               [U,S,V]=svd(subAtemp(:,:,j),0);
               avector=(diag(S));
               bvector=(avector-lambda*delta*ones(length(avector),1)).*(avector>lambda*delta*ones(length(avector),1));
               if p==q,
                  subBtemp(:,:,j)=U*diag(bvector)*V';
               else
                  BS=zeros(p,q);
                  BS(1:min(p,q),1:min(p,q))=diag(bvector);
                  subBtemp(:,:,j)=U*BS*V';
               end
           end
     Btemp=reshape(subBtemp,p,s*q);          
           %hfunction(Y,X,Btemp,lambda)
     temphvalue{t}=hmatrixfunctiondividen(totalY,totalX,B{t},lambda,s,n);
     B{t+1}=Btemp;
     alpha{t+1}=(1+sqrt(1+(2*alpha{t})^2))/2; 
     t=t+1;
end  

%estimateB1=B{t}(:,1:q);
%estimateB2=B{t}(:,(q+1):(2*q));

%profile viewer

regularizedB=B{t};
%profile viewer
%toc
end