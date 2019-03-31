function[dlossvalue]=dlossmatrixdividen(Y,X,B,s,n)
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
for j=1:s,
    for i=1:n,
        tempdloss1(:,:,j,i)=-2*X(i,j)*residual((i*row-row+1):(i*row),:);
    end
end
tempdloss=sum(tempdloss1,4);
dlossvalue=reshape(tempdloss,row,col)/n/2;
end