% ============
% Description:
% ============
%
%             This is the main program for Functional Quasi-Likelihood
%             Regression. The response is functional(densely observed) and 
%             the predictor is a vector. PACE is first used to  obtain
%             principal component scores and then SPQR is run to fit their
%             dependence on the predictor. The details are specified in
%             Chiou and M\"uller (2003). 
%             
%             
%             
%
% ======
% Usage:
% ======
%
% function [ X, W ] = FQR( y, t, z, kernel,bw,linktyp,vartyp,theta,fig, p)
% 
% Input:
%
%      y:          1*n cell array, y{i} is the vector of measurements for the ith subject,
%                  i=1,...,n.
%      t:          1*n cell array, t{i} is the vector of time points for the ith subject on which
%                  corresponding measurements y{i} are taken, i=1,...,n.
%      z:          n*p matrix, z(:,j) is the jth predictor for the ith
%                  subject
%      p:          a struct obtained from setOptions.m sets the rest of arguments for PCA.m
%                  ex:
%                  >> p = setOptions();
%
% The following input parameters are used in SPQR.
%
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
%      

% Output:
%      X:         a cell array that contains all returned values from PCA.m
%                 where the last element contains the names of X

%      W:         a K*1 cell array that contains all returned values
%                 from iterGLM.m for K principal scores




%%
function [ X, W] = FQR( y, t, z, p, fig,kernel,bw,linktyp,vartyp,theta)

if nargin<10 
    theta=[]; 
end
if nargin<9||isempty(vartyp) 
    vartyp=0; 
end
if nargin<8||isempty(linktyp) 
    linktyp=0; 
end
if nargin<7 
    bw=[]; 
end
if nargin<6||isempty(kernel) 
    kernel='epan'; 
end

if nargin<5 |isempty(fig)
    fig=1;
end
if nargin < 4
    p = setOptions();
end


X=FPCA(y,t,p);

xi_est = getVal(X,'xi_est');  

[n,K]=size(xi_est);

W=cell(K,1);
Wnames={'eta','mu','der','sigma2','beta','alpha','vbeta','pval','dis','optbw','ggrid','gfct','vgrid','vfct','names'};

for j=1:K
    [eta,mu,der,sigma2,beta,alpha,vbeta,pval,dis,optbw,ggrid,gfct,vgrid,vfct]=iterGLM(xi_est(:,j)',z,kernel,bw,linktyp,vartyp,theta,fig);
    W{j,1}={eta,mu,der,sigma2,beta,alpha,vbeta,pval,dis,optbw,ggrid,gfct,vgrid,vfct,Wnames};
end

end


