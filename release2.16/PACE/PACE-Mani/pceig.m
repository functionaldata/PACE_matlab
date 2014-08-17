% function estimating the empirical principal component basis 
%
% X,T - observed data matrix and balanced time grid 
% r - output dimension, first r pricipal components 
% xi,phi,lambda - first r pricipal components, eigenfuntions and
% eigenvalues

function [xi,phi,lambda] = pceig(X,T,r)

N = size(X,1);
M = size(X,2);
mu = mean(X);
CovX = (X'*X-N*mu'*mu)/(N-1);
CovX = (CovX+CovX')/2;   
[vec,val] = eig(CovX);   
[sval,indx] = sort(diag(val),'descend');
if nargin<3|isempty(r)
    r=1;
    % fraction of variance expained > 85%
    while sum(sval(1:r))/sum(sval)<0.85 
        r=r+1;
    end
end
if r>min(N,M) r=min(N,M); end

phi = sqrt((M-1)/range(T))*vec(:,indx(1:r))';
lambda = diag(phi*CovX*phi'*(range(T)/(M-1))^2)';
xi = (X-repmat(mu,[N,1]))*phi'/((M-1)/range(T));

end