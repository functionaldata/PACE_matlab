% 10 fold cross validation function to choose bandwidth h, using K
% nearest neighbors
%
% S - a structure, see maniMDS
% kernel - smoothing kernel
% K - selected K, as character
% Kcv - candidates and corresponding mse

function [K,Kcv] = hknn10cv(S,kernel)

Ncv = 10;

sse = zeros(1,10);
randind = randperm(S.N);
Kcan = 4:ceil((S.N/4-4)/9):floor(S.N/4);
for k=1:length(Kcan)
    for i=1:Ncv
        if i~=Ncv
            validset=randind((1:floor(S.N/Ncv))+(i-1)*floor(S.N/Ncv));
            trainset=randind([1:((i-1)*floor(S.N/Ncv)),(i*floor(S.N/Ncv)+1):S.N]);
        else
            validset=randind(((Ncv-1)*floor(S.N/Ncv)+1):S.N);
            trainset=randind(1:((Ncv-1)*floor(S.N/Ncv)));
        end
        Stmpt.Y = S.Y(trainset,:);
        Stmpt.N = length(trainset);
        Stmpt.M = S.M;
        Stmpt.X_reg = S.X_reg(trainset,:);
        Xtmpt = maniKS(S.Y(validset,:),Stmpt,kernel,num2str(Kcan(k)),0);
        sse(k) = sse(k)+trace(L2_distance(Xtmpt',S.X_reg(validset,:)',0)*sqrt(range(S.t)/(S.M-1)));
        clear Stmpt Xtmpt idxtmpt;
    end
    mse(k) = sse(k)/S.N;
end
[tmpt,id] = min(mse);
K = num2str(Kcan(id));
Kcv = {Kcan,mse};

end          
        