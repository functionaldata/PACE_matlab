function [invalid covEst] = GetRepfCov(bw,kernel,xin,yin,win,out21,out1,count)
active=find(win~=0);
xin=xin(:,active);
yin=yin(active);
win=win(active);
invalid=0;
m1 = length(out21);
m2 = length(out1);
covEst=zeros(m1,m1,m2);
for i=1:m2
  for j = 1:m1 
    for k =j:m1
     xnow = [out21(j),out21(k),out1(i)]';
     [invalid temp] = GetEstimate(xin, yin, win, bw, kernel, xnow, count); 
     if invalid == 1
          covEst = [];
          return;
     end
     covEst(j,k,i) = temp;
    end
  end
  a = triu(covEst(:,:,i),1);    %obtain upper triagular part of the covEst matrix
  covEst(:,:,i) = a+a'+diag(diag(covEst(:,:,i))); %assign the lower triangular part of the covEst matrix to be the same as the upper triangular part
end    
end