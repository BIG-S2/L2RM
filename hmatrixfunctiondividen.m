function[hfunctionvalue]=hmatrixfunctiondividen(Y,X,B,lambda,s,n)
hfunctionvalue=lossmatrixdividen(Y,X,B,s,n)+Jmatrixfunction(lambda,B,s);
end