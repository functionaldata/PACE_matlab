%function fccb = FSVDdiag(x, t_x, param_X, y, t_y, param_Y, nsvd, bwccov, alpha)
% FSVDdiag    returns 100(1-alpha)% confidence interval of functional
%             correlation, by utilizing the output of FSVD.
% Input:
% x    :  1*n cell array for function x,  x{i} is the row vector of measurements
%         for the ith subject, i=1,...,n. 
%
% t_x  :  1*n cell array, t_x{i} is the row vector of time points for the ith
%         subject at which corresponding measurements x{i} are taken,
%         i=1,...,n.
%
% y    :  1*n cell array for function y,  y{i} is the row vector of measurements
%         for the ith subject, i=1,...,n. 
% 
% t_y  :  1*n cell array, t_y{i} is the row vector of time points for the ith
%         subject at which corresponding measurements y{i} are taken,
%         i=1,...,n.
%
% param_X: an object that is returned by setOptions(), it sets the input
%         arguments for FPCA() of the X (predictor) functions (for default, set param_X = []).
%         For default values of parameters, see setOptions() for more details.
%
% param_Y: an object that is returned by setOptions(), it sets the input
%         arguments for FPCA() of the X (predictor) functions (for default, set param_X = []).
%         For default values of parameters, see setOptions() for more details.
%
% nsvd:   positive integer. Number of singular components used in regression.
%         Default nsvd = [], then use 'FVE' (fraction of variance explained) 
%         criteria to select the number of singular components that explain at least
%         75% of total variation.
%
% bwccov:  1*2 vector, bandwidths for covariance surface used for
%          smoothing of cov(X(t),Y(s))
%          bwccov(i): ith coordinate of bandwidth vector, i=1,2.
%          bwccov(1)>0 & bwccov(2)>0: user-specified bandwidths.
%          bwccov(1)==0 & bwccov(2)==0: use generalized cross-validation
%          (GCV) for automatic selection.                     [Default]
%          For the purpose of estimating singular value, the GCV-chosen 
%          bandwidth is multiplied by an empirical factor depending on the
%          sparsity of the data.
%
% alpha:  level of confidence interval, optional; 100(1-alpha)% confidence
%         interval will be created.  alpha is .05 by default.
%
% Output:
%  fccb:     a 1*2 vector, the lower bound and the upper bound of the confidence 
%            interval of functional correlation.


function [fccb fcsim] = FSVDdiag(x, t_x, param_X, y, t_y, param_Y, nsvd, bwccov, alpha)

    nsim = 200;
    
    if nargin < 9 || isempty(alpha)
        alpha = 0.05;
    end
    n = length(x);
    m = n; %max(50, floor(n/4));
    id = reshape(mysample(1:n,m*nsim,1),m,nsim);
    fcsim = zeros(1,nsim);
    tfcsim = zeros(1,nsim);
    isim = 1;
    while isim<=nsim
        xsim =  x(id(:,isim));
        txsim = t_x(id(:,isim));
        ysim = y(id(:,isim));
        tysim = t_y(id(:,isim));
        save tempsvdd xsim txsim ysim tysim fcsim;
        [res invalid]= FSVD(xsim, txsim, param_X, ysim, tysim, param_Y, nsvd, bwccov);
        if invalid == 0
            fc = getVal(res,'fc');
            fcl = 4*(fc-1/2)+1/12*(fc-1/2)^3;
            fcsim(isim) = exp(fcl)/(1+exp(fcl));
            tfcsim(isim) = getVal(res,'tfc');
            isim = isim + 1;
        end
    end
    fccb = quantile(fcsim,[alpha/2 1-alpha/2]);
    
end