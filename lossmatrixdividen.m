function[lossvalue]=lossmatrixdividen(Y,X,B,s,n)
[row col]=size(B);
subcol=col/s;

tempsum=[];
for j=1:s,
    subB{j}=B(:,((j-1)*subcol+1):(j*subcol));
    subX{j}=X(:,j);
    tempsum(:,:,j)=kron(subX{j},subB{j});
end

fit=sum(tempsum,3);
residual=Y-fit;
lossvalue=(norm(residual,'fro'))^2/n/2;
end