%===========
%Description:
%===========
%
%            This program is for generalized functional linear regression
%            with unknown link and variance functions. The response
%            variable must be a scalar and the predictor should be a
%            vector. The user can also specify a link function or a variance
%            function or both from a variety of choices provided in the
%            algorithm. The algorithm uses the iteratively reweighted least
%            squares method (IWLS) (see Chiou and M\"uller (1998), M\"uller
%            and Stadtm\"uller (2005) for details). If one deals with data
%            that include functional predictors, SPQR can be run within the
%            PACE (Principal Analysis by Conditional Expectation) program
%            where the necessary preprocessing steps are carried out. The
%            integrated PACE-SPQR program can be found at
%            http://anson.ucdavis.edu/~btliu/PACE/.
%  
%  
%======
%Usage:
%======
%
% function [eta,mu,der,sigma2,beta,alpha,vbeta,pval,dis,optbw,ggrid,gfct,vgrid,vfct] = iterGLM(Y,X,kernel,bw,linktyp,vartyp,theta,fig)
%
%======
%Input: 
%======
%      Y:          1*n vector where Y(i) is the response variable for the
%                  ith subject
%      X:          n*m matrix where X(i,:) is the predictor for the ith subject
%      kernel:     kernel function to be used for smoothing link and
%                  variance functions
%                  'epan' - Epanechnikov: 0.75*(1-x^2) on [-1,1] 
%                  'rect' - rectangle:    0.5 on [-1,1]
%                  'gauss'- gaussian:     exp(-x^2/2)/sqrt(2*pi) on [-4,4]
%                  Default is 'epan'.
%      bw:         a vector of length 3 if specified, correspondingly
%                  giving the bandwidth for smoothing the link function,
%                  the derivative of link function and the variance
%                  function 
%                  GCV will be used if not specified
%      linktyp:    an integer indicating the type of link function
%                  0 - unknown; 1 - identity; 2 - power 2; 3 - log; 4 -
%                  logit; 5 - cloglog; 6 - inverse; 7 - sqrt
%                  Default is 0.
%      vartyp:     an integer indicating the variance function
%                  0 - unknown; 1 - constant; 2 - binomial; 3 - poisson; 4
%                  - gamma
%                  Default is 0.
%      theta:      a vector of length n needed only if vartyp is 2 or 4;
%                  theta(i) givs the sample size ni for the ith subject
%                  from binomial B(ni,pi) or gamma G(ni,lambdai) distribution.
%      fig:        an integer indicating whether link and variance
%                  functions should be plotted; 0 - no, the link
%                  function and variance function evaluated on a dense grid
%                  are given in output gfct, vfct; 1 - yes
%                  Default is 1.
%
%=======
%Output:  
%=======  
%      eta:       1*n vector where eta(i) is the estimated f(E(Y|mu)) for
%                 the ith subject (f is the link function) 
%      mu:        1*n vector where mu(i) is the estimated g(eta) (g=f^{-1})
%                 for the ith subject
%      der:       1*n vector where der(i) is the estimated g'(eta) for the
%                 ith subject
%      sigma2:    1*n vector where sigma2(i) is the estimated var(Y|mu) for
%                 the ith subject
%      beta:      1*m vector of the regression coefficients; if the link
%                 function is unknown, beta is normalized to have norm 1
%                 with first component nonnegative (without intercept)
%      alpha:     estimated intercept (0 if link function is unknown) 
%      vbeta:     covariance matrix of beta (overdispersion is included,
%                 for non-overdispersed case, divide by output dis 
%                 below); if link function is known, this is the
%                 covariance matrix of [alpha,beta] 
%      pval:      1*m vector of p-values for testing beta=0; if link
%                 function is known, this is 1*(m+1) vector of p-values for
%                 testing [alpha,beta]=0
%      dis:       estimated Pearson overdispersion parameter (1 if link
%                 function is unknown)
%      optbw:     1*3 vector giving the final bandwidths for smoothing the
%                 link function, the derivative of link function and the
%                 variance function from GCV; NaN if the corresponding
%                 function is known. 
%      ggrid:     1*100 vector giving a dense time grid for gfct 
%      gfct:      1*100 vector giving the evaluated link function on time
%                 grid ggrid 
%      vgrid:     1*100 vector giving a dense time grid for vfct 
%      vfct:      1*100 vector giving the evaluated variance function on time
%                 grid vgrid   

function [eta,mu,der,sigma2,beta,alpha,vbeta,pval,dis,optbw,ggrid,gfct,vgrid,vfct] = iterGLM(Y,X,kernel,bw,linktyp,vartyp,theta,fig)

if nargin<8|isempty(fig) fig=1; end
if nargin<7 theta=[]; end
if nargin<6|isempty(vartyp) vartyp=0; end
if nargin<5|isempty(linktyp) linktyp=0; end
if nargin<4 bw=[]; end
if nargin<3|isempty(kernel) kernel='epan'; end

if (linktyp==3|linktyp==7)&length(find(Y<0))>0 error('Some responses are negative!'); end
if (linktyp==4|linktyp==5)&length(find((Y.*(1-Y))<0))>0 error('Some responses are negative or larger than 1!'); end
if linktyp==6&length(find(Y==0))>0 error('Some responses are 0!'); end
if (vartyp==2|vartyp==4)&isempty(theta) error('Theta needs to be given!'); end
if vartyp==2&length(find((Y.*(1-Y))<0))>0 error('Some responses are negative or larger than 1!'); end
if (vartyp==3|vartyp==4)&length(find(Y<0))>0 error('Some responses are negative!'); end
    
% this version is modified from the PACE integrated version
n=length(Y);
datatyp=0;
Xreg=X;
p=size(Xreg,2);
newdata=1:p;
T=[];

% initialization
if vartyp==0 | vartyp==1
    sigma2=repmat(var(Y),[1,n]);
else 
    Yadj=iniadj(Y,[],vartyp);
    sigma2=varfun(Yadj,Y,vartyp,theta);
    clear Yadj;
end
if length(bw)==3
    bw12=bw(1:2);
    bw3=bw(3);
else
    bw12=[];
    bw3=[];
end
if linktyp>0
    Yadj=iniadj(Y,linktyp,[]);
    eta=linkfun(Yadj,linktyp,0);
    der=linkfun(Yadj,linktyp,1);
    Xreg2=([ones(1,n); Xreg'])';
    clear Yadj;
    b=betaupd(Xreg2,Y,eta,Y,der,sigma2);
    eta=b*Xreg2';
    mu=linkfun(eta,linktyp,-1);
    der=linkfun(mu,linktyp,1);
    optbw1=NaN;
    optbw2=NaN;
else
    b=repmat(1/sqrt(p),[1,p]); 
    eta=b*Xreg';
    [mu,der,optbw1,optbw2]=gupd(Y,eta,sigma2,kernel,bw12);
end
[sigma2,optbw3]=sigma2upd(Y,theta,mu,vartyp,kernel,bw3);

% iteration
dbeta=1;
iter=0;
% fprintf(1,'Begin iterated reweighted least squares algorithm! \n')
while dbeta>0.01 & iter<50 
    iter=iter+1;
    b_old=b;
    if linktyp>0
        b=betaupd(Xreg2,Y,eta,mu,der,sigma2); 
        eta=b*Xreg2';
        mu=linkfun(eta,linktyp,-1);
        der=linkfun(mu,linktyp,1);
    else
        b=betaupd(Xreg,Y,eta,mu,der,sigma2);
        if datatyp==0 & length(T)>0
            bnorm=sqrt(trapz(T,b.^2));
            b=b/bnorm;
        else
            bnorm=norm(b);
            b=b./bnorm;
        end
        eta=b*Xreg';
        [mu,der,optbw1,optbw2]=gupd(Y,eta,sigma2,kernel,bw12); 
    end 
    [sigma2,optbw3]=sigma2upd(Y,theta,mu,vartyp,kernel,bw3);
    dbeta=norm(b-b_old)/norm(b_old);
end  

% if iter<50
%     fprintf(1,['Solution found after ' num2str(iter) ' iterations ! \n'])
% else
%     fprintf(1,'Maximun number of iterations (50) have been completed. Solution may not converge ! \n')
% end

if b(1)<0 b=-b;eta=-eta; end
if linktyp>0  
    if datatyp==0
        beta=b(2:end);
    else
        beta=b(2:end)*Phi;
    end
    alpha=b(1);
else
    if datatyp==0      
        beta=b;
    else
        beta=b*Phi;
    end
    alpha=0;
end

if vartyp>0
    if linktyp>0
        dis=sum((Y-mu).^2./sigma2)/(n-p-1);
    else
        dis=sum((Y-mu).^2./sigma2)/(n-p);
    end
else
    dis=1;
end

if vartyp>0
    Vinv=diag((dis*sigma2).^-1);
else
    Vinv=diag(sigma2.^-1);
end
if linktyp>0
    Dx=repmat(der',[1,p+1]).*Xreg2;
    Df=eye(p+1);  
else
    Dx=repmat(der',[1,p]).*Xreg;
    Df=(eye(p)-b'*b)/bnorm;
end
Sigmainv=n*pinv(Dx'*Vinv*Dx);
if datatyp==0
    vbeta=Df*Sigmainv*Df'/n;
elseif linktyp>0
    Phi2=zeros(p+1,size(Phi,2)+1); 
    Phi2(1,1)=1;
    Phi2(2:end,2:end)=Phi;
    vbeta=Phi2'*Df*Sigmainv*Df'*Phi2/n;
else
    vbeta=Phi'*Df*Sigmainv*Df'*Phi/n;
end

if linktyp==0
    Zsco=beta./sqrt(diag(vbeta))';
else
    Zsco=[alpha,beta]./sqrt(diag(vbeta))';
end
pval=2*(1-normcdf(abs(Zsco)));

optbw=[optbw1,optbw2,optbw3];

ggrid=min(eta):range(eta)/99:max(eta);
vgrid=min(mu):range(mu)/99:max(mu);
if linktyp>0
    gfct=linkfun(ggrid,linktyp,-1);
else
    [invalid,gfct]=locpoly(optbw(1),kernel,[],1,0,eta,mu',1./sigma2,ggrid);
    if invalid==1
        tmptbw=minbwd([eta,ggrid],3);
        clear gfct invalid;
        [invalid,gfct]=locpoly(tmptbw,kernel,[],1,0,eta,mu',1./sigma2,ggrid);
    end
end
if vartyp>1
    vfct=varfun(vgrid,[],vartyp,ones(1,100));
elseif vartyp==1
    vfct=repmat(var(Y-mu),[1,100]);
else
    [invalid2,vfct]=locpoly(optbw(3),kernel,[],1,0,mu,sigma2',ones(1,n),vgrid);
    if invalid2==1
        tmptbw=minbwd([mu,vgrid],3);
        clear vfct invalid2;
        [invalid2,vfct]=locpoly(tmptbw,kernel,[],1,0,mu,sigma2',ones(1,n),vgrid);
    end
end
if fig==1
    figure
    plot(eta,Y,'.r')
    hold on
    plot(ggrid,gfct,'-b')
    hold off
    xlabel('estimated eta')
    ylabel('observed Y')
    if linktyp==0 
        title('estimated link function') 
    else
        title('link function') 
    end
    figure
    plot(mu,(Y-mu).^2,'.r')
    hold on
    plot(vgrid,dis*vfct,'-b')
    hold off
    xlabel('estimated mu')
    ylabel('estimated squared residual')
    if vartyp==0 
        title('estimated variance function') 
    else
        title('variance function with overdispersion') 
    end
end

end
    