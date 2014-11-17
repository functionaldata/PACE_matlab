%Choose the best number of principal components through BIC method
function [no_opt,bic]=no_BIC(y,t,mu,bw_xcov,ngrid,regular,maxk,out1, out21,kernel,error,cut,rcov,xcov,npoly)
if nargin < 15
  npoly = 1;
end
bic=Inf*ones(1,maxk);
N = length(cell2mat(t));

if error == 1
   [invalid, sigma]=pc_covE(t,out1,bw_xcov,ngrid,cut,kernel,rcov, npoly);
else
   sigma = 0;
end

[lambda, phi,eigen, noeig]= getEigens(xcov,out1,out21,maxk);
[muSub, phiSub] = convertMuPhi(t, out1, mu, phi,regular);
if noeig < maxk
fprintf(1, ['Warning: max number of PC for BIC selection is no more than ' num2str(noeig) ' ! You may want to increase maxk = ' num2str(maxk) ' to a larger number to include greater flexibility.\n']);
  maxk = noeig;
end

for k=1:maxk
    sig = lambda(1:k);
    %phii = phi(:,1:k);
    
    % This part modified by Cong (Start) 10/05/2012
    [invalid, logLik] = getLogLik1(y,sig,muSub,phiSub,sigma,regular);
    if invalid
        fprintf(1, 'Warning: the covariance matrix of the estimated function is nearly singular! Reset to FVE method now!\n');
        no_opt = []; bic = [];
        return;
    end
    
    % This part modified by Cong (End) 10/05/2012
    
    bic(k)=logLik+log(N)*k; 
    tmp = bic(~isinf(bic));
    if length(tmp) > 1 &&  k > 1
       if bic(k) > bic(k-1)
           break;
       end
    end
end
bic = bic(~isinf(bic));
[temp,no_opt]=min(bic);
    
end
