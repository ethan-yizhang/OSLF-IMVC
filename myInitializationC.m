function [C] = myInitializationC(KH,k)

numker = size(KH,3);
Sigma0 = ones(numker,1)/numker;
avgKer  = mycombFun(KH,Sigma0);
[H_normalized1] = mykernelkmeans(avgKer, k);
H_normalized1 = H_normalized1./ repmat(sqrt(sum(H_normalized1.^2, 2)), 1,k);
[IDX, C] = kmeans(H_normalized1,k  , 'MaxIter',200, 'Replicates',30);
C = orth(C);