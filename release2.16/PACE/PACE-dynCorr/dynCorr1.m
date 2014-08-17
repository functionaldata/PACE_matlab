%% Dynamic corr regular case (no smoothing)

function [dyncorr corr meanmatrix]=dynCorr1(y)
%tic
% number of subjects
nsub=size(y,2);
% number of functions
ndep=size(y,1);
% number of timepoints
T=size(y{1,1},2);

% mean across subjects
meanmatrix=cell(1,ndep);
for i=1:ndep
     sub_mean=zeros(1,T);
     for k=1:T
         for j=1:nsub
             sub_mean(k)=sub_mean(k)+y{i,j}(k);
         end
     end
     sub_mean=sub_mean./nsub;
     meanmatrix{i}=sub_mean;
end

fstar=cell(ndep,nsub);
fhat=cell(ndep,nsub);
mhat=zeros(ndep,nsub);
sd=zeros(ndep,nsub);

for i=1:ndep
    for j=1:nsub
        fhat{i,j}=y{i,j}-meanmatrix{i};
        mhat(i,j)=sum(fhat{i,j})/T;
        sd(i,j)=sqrt(sum((fhat{i,j}-mhat(i,j)).^2)/T);
        fstar{i,j}=(fhat{i,j}-mhat(i,j))./sd(i,j);
    end
    i
end

% correlation
dyncorr=zeros(ndep,ndep);
corr=cell(1,nsub);
for k=1:nsub
    temp=zeros(ndep,ndep);
    for i=1:ndep
        for j=i:ndep
            temp(i,j)=sum(fstar{i,k}.*fstar{j,k})/T;
        end
        i
    end
    corr{k}=temp;
    dyncorr=dyncorr+temp;
end
dyncorr=dyncorr./nsub;
%toc
