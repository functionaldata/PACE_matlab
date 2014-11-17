% function using 10 fold cross validation to choose K, h and delta
%
% X,T - observed data cell
% datatype -integer, see maniMDS
% X_cv - data on pooled time grid, created in maniMDS 
% t -  pooled time grid
% kernel - smoothing kernel, see maniMDS for options
% Dis - distance matrix
% method - string giving the method name, see maniMDS
% d - output dimension
% Kcan, deltaca, hcan - candidate sets of the parameters if user want to
% specify them
% K,delta,h - selected parameter
% mse - mean squared errors for all candidates

function [K,delta,h,mse] = par10cv(X,T,datatype,X_cv,t,kernel,Dis,method,d,Kcan,deltacan,hcan)

if nargin<11 hcan=[]; end
if nargin<10 deltacan=[]; end
if nargin<9 Kcan=[]; end

Ncv = length(X);
if isempty(hcan)
    [Sort_dis,Ord] = sort(Dis);
    hmin = median(Sort_dis(2,:));
    hmax = max(Sort_dis(floor(size(Dis,1)*0.75)+1,:));
    hcan = hmin*(hmax/hmin).^(0:1/9:1);
end
if isempty(deltacan)&strcmp(method,'p-isomap')
    deltacan = [0,0.02,0.05,0.1,0.2];
elseif strcmp(method,'isomap')
    deltacan = 0;
end
if isempty(Kcan) Kcan = [3,5,8,12,16]; end
mse = repmat(inf,[length(Kcan),length(deltacan),length(hcan)]);
for k1=1:length(Kcan)
    for k2=1:length(deltacan)
        [tmptY,tmptIndex] = maniMethods(Dis,method,d,Kcan(k1),deltacan(k2));
        N = length(tmptIndex);
        Ntotal = 0;
        for j=tmptIndex
            Ntotal = Ntotal+length(T{j});
        end
        randind = randperm(N);
        for k=1:10
            counts = 0;
            sse = 0;
            for i=1:Ncv
                if i~=Ncv
                    validset=randind((1:floor(N/Ncv))+(i-1)*floor(N/Ncv));
                    trainset=randind([1:((i-1)*floor(N/Ncv)),(i*floor(N/Ncv)+1):N]);
                else
                    validset=randind(((Ncv-1)*floor(N/Ncv)+1):N);
                    trainset=randind(1:((Ncv-1)*floor(N/Ncv)));
                end
                Stmpt.Y = tmptY(trainset,:);
                Stmpt.N = length(trainset);
                Stmpt.M = length(t);
                Stmpt.X_reg = X_cv(trainset,:);
                Xtmpt = maniKS(tmptY(validset,:),Stmpt,kernel,hcan(k),0);
                idxtmpt = find(~isnan(Xtmpt(:,1)));
                for j=1:length(idxtmpt)
                    counts = counts+length(T{validset(idxtmpt(j))});
                    [tmpt1,tmpt2,tmptindext] = intersect(T{validset(idxtmpt(j))},t);
                    sse = sse+sum((Xtmpt(idxtmpt(j),tmptindext)-X{validset(idxtmpt(j))}).^2);
                    clear tmpt1 tmpt2 tmptindext;
                end
                clear Stmpt Xtmpt idxtmpt;
            end
            if counts/Ntotal>0.9 mse(k1,k2,k) = sse/counts; end
            % if k>1&mse(k1,k2,k)>1.1*mse(k1,k2,k-1) break; end            
        end
        clear tmptY tmptIndex randind;
    end
end
tmptmin = min(min(min(mse)));
indexnum = find(mse==tmptmin,1);
k = ceil(indexnum/(size(mse,1)*size(mse,2)));
k2 = ceil((indexnum-(k-1)*size(mse,1)*size(mse,2))/size(mse,1));
k1 = indexnum-(k-1)*size(mse,1)*size(mse,2)-(k2-1)*size(mse,1);
h = hcan(k);
delta = deltacan(k2);
K = Kcan(k1);

end
