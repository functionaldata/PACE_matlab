% ============
% Description:
% ============
%
%            This program can be used for estimating the stickiness coefficient  
%            for sparsely observed longitudinal data.
%
%            
%            The population definition is given by
%            S_X=\frac{E[\{X(T_1)-\mu(T_1)\}\{X(T_2)-\mu(T_2)\}]}
%            {[\text{Var}\{X(T_1)-\mu(T_1)\}]^\frac{1}{2}[\text{Var}\{X(T_2)-\mu(T_2)\}]^\frac{1}{2}}
%            where $\mu(t)=E[X(t)]$ is the mean trajectory and $T_1$ and $T_2$ 
%            are independent random times, independent of the process $X$, 
%            that are uniformly distributed on $\mt.$ 
%            
%
%            Reference: 
%            Gottlieb, A. and M\"{u}ller, H.G. (2012). A Stickiness Coefficient 
%            for Longitudinal Data. Computational Statistics and Data
%            Analysis (Accepted).
%            
%
% ======
% Usage:
% ======
%
%      function [stickiness_coefficient] = stick(y, t, param_X)
%
% ================
% Input Arguments:
% ================
%
%      Input y:          1*n cell array, y{i} is the vector of measurements for the
%                        ith subject, i=1,...,n. See PCA() for more details
%
%      Input t:          1*n cell array, t{i} is the vector of time points for the
%                        ith subject for which corresponding measurements y{i} are
%                        available, i=1,...,n. See PCA() for more details
%
%      param_X:          an object that is returned by setOptions() that sets the input
%                        arguments for FPCA() for X. (for default, set
%                        param_X = []). In paper param_X =
%                        setOptions('selection_k', 'FVE','FVE_threshold', 0.9);
%
% =================
% Output Arguments:  
% =================  
%      stickiness_coefficient: estimated stickiness coefficent for X
%                       
%

function [stickiness_coefficient] = stick(y, t, param_X)

p = path;
isExist = regexp(p, 'PACE');
if isempty(isExist) == 1
  addpath(genpath('../PACE/'));
end
%Run FPCA 
if (nargin == 2 |isempty(param_X))
    param_X = setOptions();
    [yy]=FPCA(y,t);
   else 
[yy] = FPCA(y,t,param_X);
end
%Extract relevant values
time_vect= getVal(yy,'out1');           %vector of time points for mu, phi and ypred
T=max(time_vect)-min(time_vect);         %gives length of time domain
lambda=getVal(yy, 'lambda');            %estimated eigenvalues
phi = getVal(yy,'phi');                 %estimated eigenfunctions



integral_vect_1=zeros(1, length(lambda));      %initialize vector to hold values of integrals
integral_vect_2=zeros(1, length(lambda));      %initialize vector to hold values of integrals

for i=1:length(lambda)
    integral_vect_1(i)=trapz(time_vect, phi(:, i));
end

for i=1:length(lambda)
    integral_vect_2(i)=trapz(time_vect, phi(:, i));
end

output=(sum(lambda.*integral_vect_1.*integral_vect_2));

stickiness_coefficient=(1/T)*output/sum(lambda);
end

