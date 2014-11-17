% ============
% Description:
% ============
% 
%   This is the main function to model volatility trajectories for high frequency observations
%   in financial markets. 
% 
% Includes the following steps:
% 
% 1) Computing scaled log-returns Z from price observations X.
% 
% 2) Computing Y process from Z, Y(t) = log(Z^2) - q_0.
% 
% 3) Functional principal component analysis for Y process to obatin functional volatility process V(t),
%    defined as V(t) = log(\sigma(t)^2).
%
% See reference:  
% M\"uller, H.G., Sen, R., Stadtm\"uller, U. (2011). Functional Data Analysis for Volatility. 
% M\"uller, H.G., Stadtm\"uller, U., Yao, F. (2006). Functional variance processes. Journal of the Americal Association, 101, 1007-1018. 
% ========
% Usage:
% ========
% 
%     [V,t, yy] = Fvola(X, T, param)
% 
%======
%Input: 
%======
%      X:          n*p matrix, X(i, j) is the price for subject i at time
%                  t(j). X must be positive and in a dense and regular
%                  design. 
%      T:          1*p vector contains t(j), 1<=j<=p. It must be a regular design, 
%                  i.e., delta= t(j+1) - t(j)are the same for all j.  
%      param:          a struct obtained from setOptions.m sets the rest of arguments for PCA.m
%                  ex:
%                  >> param = setOptions();
%=======
%Output:  
%=======  
%     V(t):        The esitmated functional volatility process V(t) = log(\sigma(t)^2).
%     t:           the coresponding time grid for V(t);
%     yy:          a cell array that contains all returned values from
%                  PCA.m for the process Y(t) = V(t) + U(t); 
%      
%
%   To see the names for yy, type names(yy)
%   To get individual value back, type getVal(yy,varname)
%   To see an example, check with example.m
%   See also FPCA, PCA, example, names
%
% Note: Squared Z are adjusted by adding the 10th percentile of % the pooled sample before taking log.
%
function [V t yy] = Fvola(X,T,param) 
    V= [];
    t = [];
    yy = [];
    n = size(X, 1);
    p = size(X, 2);
    delta = T(2:p)- T(1:(p-1));
    if (min(min(X))<= 0 || (min(delta) ~= max(delta)))
     fprintf(1,'Error: the input value of X must be positive and sampled at equally spaced time points.');
     return;
    end
    delta = delta(1);
    Z  = 1/sqrt(delta)*log(X(:, 2:p)./X(:, 1:(p-1)));
    threshold2 =  quantile(reshape(Z.^2, n*(p-1), 1), 0.1);
    Y = log(Z.^2 + threshold2) + 1.27; % here 1.27 is a constant q_0 defined in the paper. 
    y = mat2cell(Y, ones(n, 1), (p-1))';
    t = mat2cell(repmat(T(1:(p-1)), n, 1), ones(n, 1), (p-1))';
    yy = FPCA(y,t,param);   
    V = getVal(yy,'y_pred');
    t = getVal(yy, 't');
end
