function beta = updateBeta_OSLFIMVC(HP,WP,Y,C)
numker = size(WP,3);
HHPWP = zeros(numker,1);
Hstar = Y*C;
for  p=1:numker
    HHPWP(p) = trace(Hstar'*(HP(:,:,p)*WP(:,:,p)));
end
beta = HHPWP./norm(HHPWP);
beta((beta<eps))=0;
beta = beta./norm(beta);