function [fhat,xihat,noeig,loglike,fhist] = getfhat(x,xrange,f_mu,phihat,xout,k,noeigtmp,kern,nbins)

    nn = length(x);
    tt = zeros(1,nn);
    for i = 1:nn
        tt(i) = length(x{i});
    end
    if k >= 1
        
        if k > noeigtmp;
            k = noeigtmp;
            fprintf(1,['Warning: at most ' num2str(noeigtmp) ' number of PC can be selected!\n']);
        end
            
        % Estimate the k FPC scores xi
        phitmp1 = cell(1,nn);
        phitmp2 = zeros(k,nn);   % k*nn
        for i = 1:nn
            phitmp1{i} = interp1(xout',phihat(:,1:k),x{i}','spline')';
            phitmp2(:,i) = mean(phitmp1{i},2);               
        end
        clear phitmp1;
        
        phi_mean = trapz(xout,repmat(f_mu',1,k).*phihat(:,1:k),1);  % 1*k

        xihat = phitmp2'-repmat(phi_mean,nn,1);     % nn*k
        fhattmp = repmat(f_mu,nn,1)'+phihat(:,1:k)*xihat';   % xout*nn
        % Remove the negative values (if any) in the tmp density
        % Adjust the tmp density to get a real density, i.e., integrates to 1
        fhattmp(fhattmp < 0) = 0;
        fhat = fhattmp*diag(1./trapz(xout,fhattmp,1));

        noeig = k;
        loglike = [];
        fhist = [];
                     
    else
        
        phitmp1 = cell(1,nn);
        phitmp2 = zeros(noeigtmp,nn);   % noeigtmp*nn
        for i = 1:nn
            phitmp1{i} = interp1(xout',phihat,x{i}','spline')';
            phitmp2(:,i) = mean(phitmp1{i},2);               
        end
        
        phi_mean = trapz(xout,repmat(f_mu',1,noeigtmp).*phihat(:,1:noeigtmp),1);   % 1*noeigtmp
        xihataic = phitmp2-repmat(phi_mean,nn,1)';      % noeigtmp*nn
        
        fhataic = cell(1,noeigtmp); 
        loglike1 = zeros(1,noeigtmp);    

        % Calculate the histogram densities
        hgrid = min(xrange):range(xrange)/(nbins-1):max(xrange);
        fhist = zeros(nn,nbins);
        for i = 1:nn
            [fhist(i,:)] = hist2fx(x{i},nbins,[],0,kern,hgrid,xrange);
        end        
        fhist = fhist';
        
        fmu = interp1(xout,f_mu,hgrid,'spline');
        phi = interp1(xout',phihat,hgrid','spline');
        for j = 1:noeigtmp
            
            fhatmp = repmat(fmu,nn,1)'+phi(:,1:j)*xihataic(1:j,:);   % nbins*nn
            fhatmp(fhatmp < 0) = 0;
            fhataic{j} = fhatmp*diag(1./trapz(hgrid,fhatmp,1));
            loglike1(j) = fdev(fhataic{j},fhist);
            
        end
        
        loglike = loglike1+2*(1:noeigtmp);
        noeig = find(loglike==min(loglike));
        
        xihat = xihataic(1:noeig,:)';
        fhat = interp1(hgrid,fhataic{noeig},xout,'spline');
        fhat(fhat < 0) = 0;
        fhat = fhat*diag(1./trapz(xout,fhat,1));
        
        fhist = interp1(hgrid,fhist,xout,'spline');
        fhist(fhist<0) = 0;
        fhist = fhist*diag(1./trapz(xout,fhist,1));
        
    end
