% This is the function to obtain the local linear smoothing estimate for
% one point xnow.

function [invalid mu] = GetEstimate(xin, yin, win, bw, kernel, xnow, count)
    %count: the number of repeats for each unique point in xin.
   
    %locating local window
    invalid = 0;
    m0 = size(xin, 2);
    d = size(xin, 1);
    if strcmp(kernel,'gauss') ~= 1      %when it is not gaussian kernel, look for the grids that have domain -1 and 1
        list = sum(abs(xin- (repmat(xnow, 1, m0))) <= repmat(bw', 1, m0) + 10^(-6), 1);
        ind = find(list ==d);    
    else
        ind = 1:m0;
    end
    clear list;
    lx=xin(:,ind);
    ly = yin(ind);
    lw=win(ind);
    if length(unique(lx','rows'))>=d
       m = size(lx, 2);
       indd=(1:m)';
       lwmat = repmat(lw, d, 1);
       llx = (lx - repmat(xnow, 1, m))./repmat(bw', 1, m);  
       if strcmp(kernel,'epan')
          temp=prod(lwmat.*(1-llx.^2)*(3/4), 1);
       elseif strcmp(kernel,'rect')
          temp=prod(lwmat.*repmat(1,d,m)/2, 1);
       elseif strcmp(kernel,'gauss')
          temp = prod(lwmat.*(exp(-0.5*(llx.^2))/sqrt(2*pi)), 1);
       elseif strcmp(kernel,'gausvar')
          temp = prod(lwmat.*(exp(-0.5*(llx.^2))/sqrt(2*pi).*(1.25-0.25*llx.^2)), 1);
       elseif strcmp(kernel,'quar')
          temp = prod(lwmat.*((1-llx.^2).^2).*(15/16), 1);
       end
       temp = temp.*count(ind);
       W=sparse(indd,indd,temp);
       % computing design matrix
       X=zeros(m,d+1);
       X(:,1)=ones(m,1);
       X(:,2:(d+1))=lx'-repmat(xnow', m, 1);
       beta=pinv(X'*W*X)*X'*W*ly;
       clear X W;
       mu = beta(1);         
    else
       invalid=1;
       fprintf(1,'No enough points in local window, please increase bandwidth\n');
       mu = [];
       return;
    end
end

