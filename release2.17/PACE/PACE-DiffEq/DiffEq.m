% ============
% Description:
% ============
%
%            This program can be used for estimating an empirical first order stochastic 
%            ordinary differential equation with time-varying coefficents
%            and a smooth drift process. It requires a single input object
%            XXd which is the output from FPCAder.
%            
%            The model is given by
%            X^(1)(t)-\mu^(1)(t)=\beta(t)[X(t)-\mu(t)]+Z(t) 
%            where \beta(t) is a time varying linear coefficent and Z(t) represents a
%            random drift process.
%
%            Reference: 
%            M\"{u}ller, H.G., Yao, F. (2010). Emprical dynamics for longitudinal data. 
%            Annals of Statistics 38, 3458-3486.
%
% ======
% Usage:
% ======
%
%      function [toutx beta1 zdc Gz R2 phiz VarXder VarZ] = DiffEq(XXd)
%
% ================
% Input Arguments:
% ================
%
%      Input XXd: the returned values from FPCder. See FPCder() for more details
%
% =================
% Output Arguments:  
% =================  
%      toutx:     1*N vector, distinct input time points with ascending
%                 order from all subjects used in XXd.   
%
%      beta1:     1*N vector, estimated dynamic tranfer function \beta(t)
%                 on grid toutx.       
%
%      zdc:       1*N cell, where zdm{i} is the estimated  drift function 
%                 evaulated on grid toutx for the ith subject. 
%
%      GZ:        N*N matrix, smoothed covariance surface for the drift
%                 process Z(t) on grid toutx
%
%      R2:       1*N vector, estimated coefficent of determination
%                 function R^2(t) on grid toutx. See Reference for more details. 
%
%      phiz:     N*K matrix, estimated principal component functions for
%                 the drift process Z(t) on grid toutx where K is the
%                 number of included eigen-components in XXd.
%
%      VarXder:   N*1 vector, estimated variance function for X^1(t) on
%                 grid toutx
%
%      VarZ:      N*1 matrix, estimated variance function for Z(t) on grid
%                 on toutx

function [ toutx beta1 zdm Gz R2 phiz VarXder VarZ] = DiffEq( XXd )




%Extract outputs from FPCAder
no_optx=getVal(XXd,'no_opt');
out21x=getVal(XXd,'out21');
toutx=getVal(XXd,'out1');
lamx=getVal(XXd,'lambda');
mufull=getVal(XXd, 'mu');
    mux=mufull{1};
    mud1=mufull{2};
phifull=getVal(XXd,'phi');
    phix=phifull{1};
    phid1=phifull{2};
x1full=getVal(XXd,'y_pred'); 
    x1=x1full{1};
    xd1=x1full{2};
    n=length(x1);
   xcovd=getVal(XXd,'xcov');
xcov0=xcovd{1};
xcovd1=xcovd{2};

%Dynamic Transfer Function
beta1=(lamx*(phid1.*phix)')./(lamx*(phix.*phix)'); 

xdd={}; zd={}; xz0={}; xz1={};
x1m=[]; xd1m=[]; zdm=[]; xz0m=[]; xz1m=[]; zm=[]; xz2m=[];
for i=1:n
    xdd{i}=mud1+beta1.*(x1{i}-mux); 
    zd{i}=xd1{i}-xdd{i}; 
    x1m(i,:)=x1{i}; xd1m(i,:)=xd1{i}; zdm(i,:)=zd{i}; 
end

%centered subject specific drift procecess 
zdmean=repmat(mean(zdm),n,1);
zdmc=zdm-zdmean; 
tz={}; zdc={};
for i=1:n
    tz{i}=toutx; 
    zdc{i}=zdmc(i,:); 
end

%Covariance Function for drift process
G00=phix*diag(lamx)*phix'; 
G11=phid1*diag(lamx)*phid1'; 
G01=phix*diag(lamx)*phid1'; G10=G01'; 
Gz_cross=(repmat(beta1',1,no_optx).*phix)*diag(lamx)*phid1'; 
Gz_last=(repmat(beta1',1,no_optx).*phix)*diag(lamx)*(phix.*repmat(beta1',1,no_optx))'; 
Gz=G11-Gz_cross-Gz_cross'+Gz_last;

%Coefficent of Determination
R2=1-diag(Gz)./diag(G11); 

%Eigen-decomposition of drift process
zcov0=zdmc'*zdmc; varz0=diag(zcov0)'; 
zcov0=(zcov0+zcov0')./2; 
Gz=(Gz+Gz')./2;
indp=find(out21x>=0);
[az,bz]=meshgrid(out21x(indp),out21x(indp));
[az0,bz0]=meshgrid(toutx,toutx);
Gzp=interp2(az0,bz0,Gz,az,bz,'spline');
Gzp=(Gzp+Gzp')./2;
[lamz, phiz, temp1, temp2] = getEigens(Gzp,toutx,out21x,no_optx,1); 

%Variance Functions of X^(1)(t) and Z(t)
VarXder=diag(G11);
VarZ=diag(Gz);

end

