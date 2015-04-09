% Function fdev calculates the deviance variance.

function [loglike] = fdev(fnew, fhist)

    if iscell(fnew) == 1
    
      nsub = length(fnew);
      tmp = zeros(1,nsub);
      for kk = 1:nsub
          tmp(kk) = functmp(fnew{kk}, fhist{kk});
      end
      loglike = 2*sum(tmp);
    
    else
    
      nsub = size(fnew,2);
      tmp = zeros(1,nsub);
      for kk = 1:nsub
          tmp(kk) = functmp(fnew(:,kk), fhist(:,kk));
      end
      loglike = 2*sum(tmp);
      
    end
 