clear;
clc;
path = '*';
addpath(genpath(path));
dataName = 'wdbc'
load(['datasets\',dataName,'_Kmatrix'],'KH','Y');
numclass = length(unique(Y));
Y(Y<1) = numclass;
numker = size(KH,3);
num = size(KH,1);
KH = kcenter(KH);
KH = knorm(KH);
epsionset = [0.05:0.05:0.5];
respath = [path];
for ie =1:length(epsionset)
    for iter = 1 : 10
        fprintf('%s: missing_ratio: %d ,iter: %d\n',dataName,epsionset(ie)*100,iter);
        load([path,'incompleteKernelDatasets\',dataName,'\',dataName,'_missingRatio_',num2str(epsionset(ie)),...
            '_iter_',num2str(iter),'.mat'],'S');
        qnorm = 2;
 
        iseedset = [0:19];
        lambdaset9 = 2.^[-15:1:15];
        res_allmean9 = zeros(4,length(lambdaset9),length(iseedset));
        res_allstd9 = zeros(4,length(lambdaset9),length(iseedset));
        Sigma_all9 = zeros(numker,length(lambdaset9),length(iseedset));
        for is = 1:length(iseedset)
        s=RandStream('mt19937ar','Seed',iseedset(is));
        RandStream.setGlobalStream(s);
         tic
            for il=1:length(lambdaset9)
            [H_normalized9,C9,WP9,Sigma_all9(:,il,is),obj9] = OS_LF_IMVC_alg(KH,S,numclass,lambdaset9(il));
            [res_allmean9(:,il,is),res_allstd9(:,il,is)] = myNMIACCV2(H_normalized9,Y,numclass);
            end
        [~,max_idx]=max(mean(res_allmean9(:,:,is)),[],'all','linear');
        res_mean9(is,:) = res_allmean9(:,max_idx,is);
        Sigma99(is,:) = Sigma_all9(:,max_idx,is);
        timecost9(is) = toc;
        end
        res_mean(:,9) = mean(res_mean9);
        res_std(:,9) = std(res_mean9);
        timecost(9) = mean(timecost9);

       
    end
end
