function[Jvalue]=Jmatrixfunction(lambda,B,s)
[row col]=size(B);
subcol=col/s;

for j=1:s,
    subB{j}=B(:,((j-1)*subcol+1):(j*subcol));
    [U,S,V]=svd(subB{j},0);
    Jvalue1(j)=sum(diag(S));
end    
    
Jvalue=lambda*sum(Jvalue1);
end