clear;

% sample size
n=100;

% the number of random permutations (decoupling) used for screening threshold selection
perrep=10;

% image dimensions
p=64;
q=64;

% total dimension of X
s=2000;

% Number of nonzero components
s0=4;

% error correlation between voxels
arnumber=0.5;

% error variance
errorvariance=1;

trueindex=[1 2 3 4];

% true coefficient matrices
matrixB1=zeros(p,q);
 matrixB1(31:34,17:30)=1;
 matrixB1(17:48,31:34)=1;
 matrixB1(31:34,35:48)=1;

matrixB2=zeros(p,q);
matrixB2(17:48,17:48)=1;

load('matrixB3.mat');

load('matrixB4.mat');

truecoefficientB=[];

truecoefficientB(1,:,:)=matrixB1;
truecoefficientB(2,:,:)=matrixB2;
truecoefficientB(3,:,:)=matrixB3;
truecoefficientB(4,:,:)=matrixB4;
truecoefficientB(5:s,:,:)=0;

% generate correlation matrix among X
sigma=zeros(s,s);
for i=1:s,
    for j=1:s,
        sigma(i,j)=0.5^(abs(i-j));
    end
end

RR = chol(sigma);


% generate correlation matrix between pixels of the random error 
pq=p*q;
errorsigma=zeros(pq,pq);

for i=1:pq,
    for j=1:pq,
        integeri=floor((i-1)/p)+1;
        remainderi=i-p*(integeri-1);
        integerj=floor((j-1)/p)+1;
        remainderj=j-p*(integerj-1);
        errorsigma(i,j)=arnumber^(abs(integeri-integerj)+abs(remainderi-remainderj));
    end    
end

sigmaRR=chol(errorsigma);

%%% generate data
RandStream.setGlobalStream (RandStream('mt19937ar','seed',2018));

Y=zeros(n*p,q);

for i=1:n,
    XX = randn(1,s);
    Xvector= XX*RR;
    X(i,:)=Xvector;
    tempX(i,:)=Xvector;
    for j=1:s0,
        tempcomponent(:,:,j)=Xvector(j)*truecoefficientB(j,:,:);
    end
    fitsum=sum(tempcomponent,3);
    error = reshape(reshape(errorvariance*randn(p,q),1,pq)*sigmaRR,p,q);
    Y((i*p-p+1):(i*p),:)=fitsum+error;
    tempY((i*p-p+1):(i*p),:)=Y((i*p-p+1):(i*p),:);
    arrayY(i,:,:)=Y((i*p-p+1):(i*p),:);
    sumtempY(:,:,i)=Y((i*p-p+1):(i*p),:);
end  

%% selceted index of screening step
includeindex=screeningmatrix(X,arrayY,p,q,s,n,perrep);

%% random seed for cross validation
RandStream.setGlobalStream (RandStream('mt19937ar','seed',2018+1000));

%% generate the random split for 5-fold cross validation
k=5;
cvsplit = cvpartition(n,'kfold',k);

%% Center the response and standardize the covariates
   meantempY=mean(sumtempY,3);
    
    meantempX=mean(tempX,1);   
    
    stdtempX=std(tempX,1);
    
    for i=1:n,
        centerY((p*i-p+1):(p*i),:)=tempY((p*i-p+1):(p*i),:)-meantempY;
        standardizedX(i,:)=(tempX(i,:)-meantempX)./stdtempX;
    end
 
  %%% selected centered X after screening
    selectstandardizedX=standardizedX(:,includeindex);

  %%% the final estimate of B corresponding to the selected X 
    regularizedB=regularizedmatrixcvestimation(selectstandardizedX,centerY,cvsplit,p,q,length(includeindex),n);
    
  selecteds=length(includeindex);
  
  %%% estimatedB is a list of objects with 
  estimatedB=[];
  for j=1:s,
      estimatedB{j}=zeros(p,q);
  end    
  
  for j=1:selecteds,
      estimatedB{includeindex(j)}=regularizedB(:,(1+q*(j-1)):(j*q));
  end   
  
  %%% plot the coefficient matrices
  %imagesc(estimatedB{1})
  %imagesc(estimatedB{2})
  %imagesc(estimatedB{3})
  %imagesc(estimatedB{4})
  %imagesc(estimatedB{5})
  
    
  