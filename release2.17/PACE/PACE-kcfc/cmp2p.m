function [adjustRand,correct_rate,class_order,ctable] = ...
    cmp2p(clust_result,clust_external)

% Input :
%    clust_result   -  nobs x 1 vector of cluster labels for the first 
%                      partition. (number of clusters = nclust)
%    clust_external -  nobs x 1 vector of cluster labels for the second
%                      partition. (number of clusters = nclass)
% Output 
%    adjustRand   -  adjusted Rand index 
%    correct_rate -  correct classification rate
%    class_order  -  nclust x 1 vector of the cluster orders 
%                    If the jth cluster of clust_result corresponds to the 
%                    ith cluster of clust_external, then class_order(j)=i.
%    ctable       -  nclass x nclust contingency table, the nclass rows 
%                    correspond to the partition of clust_external and the
%                    nclust columns correspond to clust_result.

nobs = length(clust_result);

[uniqueclust_r,indx_a1,indx_b1] = unique(clust_result);
[uniqueclust_e,indx_a2,indx_b2] = unique(clust_external); 
nclust = length(uniqueclust_r);
nclass = length(uniqueclust_e);
 
Amat = zeros(nobs,nclass); 
Bmat = zeros(nobs,nclust); 
for i = 1:nobs;
    Amat(i,indx_b2(i)) = 1;
    Bmat(i,indx_b1(i)) = 1;
end;
ctable = Amat'*Bmat;
rowsum = sum(ctable,2);
colsum = sum(ctable,1);
total = sum(sum(ctable));

index = 0; 
for i = 1:nclass;
    for j = 1:nclust;
        if ctable(i,j) >= 2;
           index = index + nchoosek(ctable(i,j),2);
        end;
    end;
end;
combclass = 0; 
for i = 1:nclass;
    if rowsum(i) >= 2;
       combclass = combclass + nchoosek(rowsum(i),2);
    end;
end;
combclust = 0;
for i = 1:nclust;
    if colsum(i) >= 2;
       combclust = combclust + nchoosek(colsum(i),2);
    end;
end;
index_expect = combclass*combclust/nchoosek(total,2);
index_max = (combclass + combclust)/2; 
adjustRand = (index - index_expect)/(index_max-index_expect);
  
if nclass ~= nclust;
   disp(['Warning: nclass ~= nclust, cannot calculate correct_rate']);
   correct_rate = -99;
   class_order = zeros(nclust,1);
else;
   nc = nclass;
   permutways = perms(1:nc);
   nperms = size(permutways,1);
   ncorrect = zeros(nperms,1);
   for k = 1:nperms;
        for j = 1:nc;
            ncorrect(k) = ncorrect(k) + ctable(permutways(k,j),j);
        end; 
    end; 
    correct_rate = max(ncorrect)/total;
    class_order =  permutways(find(ncorrect == max(ncorrect)),:)';
end;