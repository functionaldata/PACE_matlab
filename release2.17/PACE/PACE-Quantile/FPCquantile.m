% ============
% Description:
% ============
% 
% This is the main function to implement the conditional Quantile estimation with functional covariates, 
% where the predictor is a function X(t_x) and the response is a scalar Y. 
% 
% Reference: Chen, K., Müller, H.G. (2011). Conditional quantile analysis
% when covariates are functions, with application to growth data.
% J. Royal Statistical Society B.
% 
% 
% IMPORTANT: This model is intended primarily for densely sampled trajectories
% 
% Includes the following steps:
% 
% 1) FPCA using the PACE method for X(t_x)
% 
% 2) Computation of the conditional distribution function through a
%    functional generalized linear model.
% 
% 3) Prediction of quantiles for given predictor values
%
% 
% ========
% Usage:
% ========
% 
% [predQ, predF, xx] = FPCfam(x, t_x, y, outQ, param_X, isNewSub)
% 
% =======
% Input:
% =======
% 
% x    :   1*n cell array for predictor function x, where x{i} is the row vector of
%          measurements for the ith subject, i=1,...,n. It may contain data for subjects
%          that are used for prediction; this is controlled by "isNewSub" which is either a vector consisting of
%          0's and 1's according to whether subject is used for prediction (0) or estimation (1), or is controlled 
%          by a positive integer nn. In this case, nn is the number of subjects to be used for estimation and 
%          n-nn is the number of remaining subjects to be used for prediction, corresponding to 
%          the last n-nn data rows. When "isNewSub" is set to [], all n subjects
%          are used for estimation and no prediction will be calculated; see "isNewSub" for more details.
% 
% t_x  :   1*n cell array, t_x{i} is the row vector of time points for the ith
%          subject at which corresponding measurements x{i} are taken,
%          i=1,...,n. It contains subjects that are used for prediction.
%          See above for 2 different cases of "isNewSub" and the definition of
%          "isNewSub" for more details.
% 
% y    :   i) When no prediction is requested, that is, isNewSub = [] or
%             isNewSub = 0:
%             a 1*n vector for scalar response y, then y(i) is the response value
%             for the ith subject, i = 1,...,n.
%          ii) When prediction is requested and y is 1*n, it will be
%              truncated to 1*nn according to the isNewSub definition below.
%             (see "isNewSub" for more details).
% outQ   : a vector of desired quantile levels. If set to [], the default
%          value outQ = [0.1, 0.25, 0.5, 0.75, 0.9] will be used. 
% param_X: an object that is returned by setOptions() that sets the input
%          arguments for FPCA() of the x (predictor) functions (for default, set param_X = []).
% 
% 
%
% isNewSub: i) 1*n vector of 0s or 1s, where
% 
%              1 : the data for the corresponding subject, i.e.,
%                  X(isNewSub == 1), t_x(isNewSub==1)
%                  are used for prediction only;
% 
%                 The count is n-nn for subjects with isNewSub = 1.
% 
%              0 : the data for the corresponding subject, i.e.,
%                  X(isNewSub == 0), t_x(isNewSub == 0), y(isNewSub == 0)
%                  are used for estimation only.
% 
%                  The count is nn for subjects with isNewSub = 0
% 
%           This option is convenient for computing leave-one-out
%           prediction if desired.
% 
% 
%          ii) If it is a positive integer, say 5, then the last 5 subjects (in the order from
%              top to bottom within the array x) and their values for
%              t_x are used for prediction. In other words, when choosing
%              this option, one would append the ``new'' subjects for which one desires prediction of the response
%              to occur at the end of x, t_x. Then the first nn rows of the arrays x, t_x and y are used 
%              for estimation and the remainder (last n-nn rows) of x,t_x will be used for prediction.

%           This option is convenient for obtaining predictions for a set of new subjects.
% 
%           iii) set isNewSub = [] for the default value, which is no prediction.
% 
% Details: i) There are no default values for the first 3 arguments, that
%             is, they are always part of the input arguments.
%         ii) Any unspecified or optional arguments can be set to "[]" for
%             default values;
%          v) When isNewSub is set to be [], no prediction will be performed,
%             then predQ and predF will contain the fitted quantiles and
%             fitted conditional distribution functions for the training
%             data.
% =======
% Output:
% =======
% 
%  predQ: a matrix of n*length(outQ): the first nn rows containing fitted
%         conditional quantiles of Y corresponding to the trainning subjects, and the last n-nn rows 
%         containing predicted conditional quantiles of Y corresponding to the subjects X(isNewSub ==1).
%       
%  predF: a matrix of n*100. The ith row contains the fitted or predicted
%         conditional distribution function F(y|X_i), evaluated at an equally
%         spaced grid of 100 points.
%  b    : a matrix of K*100 contains the coefficient functions, defined as 
%         F(y|X) = g(\sum_(k=0)^K b_k(y)\xi_k), see equation (5) in the paper for details,
%         where K is the number of components selected to expand the predictor functions X, and 
%         \xi_k  is the kth principal component score.
%  xx:   an aggregated object that contains the returned values from
%        FPCA(x,t_x,param_X).
%        See PCA() or pcaHELP.txt for more details. Use getVal() to retrieve
%        the individual values.
% 
%    o    See exampleFPCquantile.m for an example of the conditional
%    quantile estimation with functional covariates.
% 

function [predQ, predF, b, xx] = FPCquantile(x, t_x, y, outQ, param_X, isNewSub)

   [x, t_x, y, t_y, isYFun, newx, new_tx, newy, new_ty, invalid] = pre(x, t_x, y, [], isNewSub);
   nn = length(x);
   n = length(x) + length(newx);
   if invalid == 1
       return;
   end

   if isempty(param_X)
       param_X = setOptions('selection_k','BIC1');   %set default for param_X
   end
       
   verbose_x = getVal(param_X,'verbose');
   if strcmp(verbose_x, 'on') == 1
       fprintf(1, 'Obtain functional object for x:\n');
   end
   
   if isempty(outQ)
       outQ = [0.1, 0.25, 0.5, 0.75, 0.9];
   end
   
   %perform PCA on x
   xx = FPCA(x, t_x, param_X);   
   xitrain=getVal(xx,'xi_est');
   
   %glm
   zm = 50;
   ysort = sort(y);
   zgrid= linspace(ysort(5), ysort(nn-5), zm);
   K_x = getVal(xx, 'no_opt');
   b=zeros(zm, K_x+1);
   for i=1:zm
       z=zgrid(i);
       ytrain =zeros(1,nn);
       ytrain(1, y <z| y ==z)=1;
       warning('off', 'stats:glmfit:IterationLimit')
       b(i,:)=glmfit(xitrain, ytrain', 'binomial','link','logit');
   end
   
   xi = zeros(n, K_x);
   xi(1:nn, :) = xitrain;
   if ~isempty(newx)
       [ypred,newpcx] = FPCApred(xx,newx,new_tx);
       clear ypred
       xi((nn+1):n, :) = newpcx; 
   end
   
   % predQ and predF 
   m = 100;
   outF = linspace(ysort(5), ysort(nn-5), m);
   predF = zeros(n, m);
   nQ = length(outQ);
   predQ = zeros(n, nQ);
   
   bw = 4*(zgrid(2)-zgrid(1));
   for i = 1:n
       F = 1:zm;
       for j =1:zm
           F(j)=glmval(b(j,:)',xi(i,:),'logit');
       end
       win=ones(1,zm);
       [invalid, predF(i, :)]=lwls(bw,'epan', 1 ,1, 0, zgrid, F', win, outF);
   end
   
   for i= 1:n
       for j = 1:nQ
           if  predF(i, m)< outQ(j) || predF(i, m) == outQ(j) 
               predQ(i, j)= outF(m);
           else
               predQ(i, j)= outF(find(predF(i, :) > outQ(j), 1, 'first'));
           end
       end
   end
end
