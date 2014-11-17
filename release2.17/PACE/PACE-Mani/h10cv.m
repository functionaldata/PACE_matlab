% 10 fold cross validation function to choose bandwidth h, using fixed
% radius
%
% S - a structure, see maniMDS
% kernel - smoothing kernel
% hmin - optional minimum bandwidth
% h - selected bandwidth
% Hcv - candidates and corresponding mse

function [h,Hcv] = h10cv(S,kernel,hmin)

if nargin<3 hmin = []; end

Ncv = 10;
[Sort_dis,Ord] = sort(S.Dis);

if isempty(hmin)
	hmin = 2*median(Sort_dis(2,:));
end

hmax = max(median(Sort_dis(floor(S.N*0.75)+1,:)),2*hmin);

hcan = hmin*(hmax/hmin).^(0:1/9:1);
sse = zeros(1,10);
counts = zeros(1,10);
randind = randperm(S.N);
for k=1:10   
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
        Xtmpt = maniKS(S.Y(validset,:),Stmpt,kernel,hcan(k),0);
        idxtmpt = find(~isnan(Xtmpt(:,1)));
        counts(k) = counts(k)+length(idxtmpt);
        sse(k) = sse(k)+trace(L2_distance(Xtmpt(idxtmpt,:)',S.X_reg(validset(idxtmpt),:)',0)*sqrt(range(S.t)/(S.M-1)));
        clear Stmpt Xtmpt idxtmpt;
    end
    if counts(k)/S.N<0.9 
        mse(k) = NaN;
    else
        mse(k) = sse(k)/counts(k);
    end
end
[tmpt,id] = min(mse);
h = hcan(id);
Hcv = {hcan,mse,counts};

end          
        