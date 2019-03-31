%%% n is the sample size, p, q are the dimensions of the 2D image, s is the dimension of the covariates
%%% perrep is the number of random permutations (decoupling) used for screening threshold selection  
%%% X is n by s, arrayY is n by p by q. 
function[includeindex]=screeningmatrix(X,arrayY,p,q,s,n,perrep)

%%% screening procedure
corrB=[];
firstsingularvalue=[];
L1normcorrB=[];
FnormcorrB=[];
vectorizeY=reshape(arrayY,n,p*q);
for l=1:s,
    corrB(l,:)=corr(X(:,l), vectorizeY); 
    S=svd(reshape(corrB(l,:),p,q));
    firstsingularvalue(l)=S(1);
end    

% random decoupling to obtain the threshold
for tt=1:perrep, 
    randomindex=randsample(n,n);
    for l=1:s,
        nullcorrB(l,:)=corr(X(:,l), vectorizeY(randomindex,:)); 
        nullS=svd(reshape(nullcorrB(l,:),p,q));
        nullfirstsingularvalue(l)=nullS(1);
    end 
    allnullvalue=nullfirstsingularvalue;  
   
    perthresholdvalue(tt)=max(allnullvalue);
end

thresholdvalue=median(perthresholdvalue);

includeindex=find(firstsingularvalue>thresholdvalue);

end
