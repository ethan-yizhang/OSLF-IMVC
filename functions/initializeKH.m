function KH = initializeKH(KH,S)

numker = size(KH,3);
num = size(KH,1);
for p = 1 : numker
    mis_set = S{p}.indx;
    obs_set = setdiff(1:num,mis_set);
    Wp = zeros(length(obs_set),length(mis_set));
    KH(obs_set,obs_set,p) = KH(obs_set,obs_set,p);
    KH(obs_set,mis_set,p) = KH(obs_set,obs_set,p)*Wp;
    KH(mis_set,obs_set,p) = KH(obs_set,mis_set,p)';
    KH(mis_set,mis_set,p) = KH(mis_set,obs_set,p)*Wp;
end