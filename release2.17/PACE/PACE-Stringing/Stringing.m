% ============
% Description:
% ============
% 
% This is the main function to implement the Stringing method. Stringing of high
% dimensional data is implemented with distance-based metric Multidimensional Scaling,
% mapping high-dimensional data to locations on a real interval, such that predictors that
% are close in a suitable sample metric also are located close to each other on the interval. 
% Established techniques from Functional Data Analysis can
% be applied for further statistical analysis once an underlying stochastic process and the
% corresponding random trajectory for each subject have been identified. 
% 
% Reference: Chen, K., Chen, K., Müller, H.G., Wang, J.L. (2011).
% Stringing high-dimensional data for functional analysis.
% J. American Statistical Association 106, 275-284.
% 
% ========
% Usage:
% ========
% 
% [newXmat, stringed, disMatrix, x, t_x, SubId, xx, pred_x, pred_pc] = Stringing(Xmat, dis_option, disMat, isNewSub, isFPCA, isPlot)
% 
% =======
% Input:
% =======
% Xmat    :n *p data matrix, where xmat(i, :) is the row vector of
%          measurements for the ith subject, i=1,...,n. It may contain data for subjects
%          that are used for prediction; this is controlled by "isNewSub" which is either
%          a vector consisting of 0's and 1's according to whether subject is used for prediction (0)
%          or estimation (1), or is controlled by a positive integer nn. In this case,
%          nn is the number of subjects to be used for estimation and n-nn is the number of
%          remaining subjects to be used for prediction, corresponding to 
%          the last n-nn data rows. When "isNewSub" is set to [], all n subjects
%          are used for estimation and no prediction will be calculated;
%          see "isNewSub" for more details.
% dis_option: one of the following:
%          'euclidean'   - Euclidean distance (default)
%          'correlation' - One minus the sample linear correlation between
%                          observations (treated as sequences of values).
%          'spearman'    - One minus the sample Spearman's rank correlation
%                          between observations (treated as sequences of values).
%          'hamming'    - Hamming distance, percentage of coordinates that differ.
%          'user'        - a user provided dissimilarity matrix named 'disMat' will be
%                          used; in this case the input 'disMat' must be provided.
% 
% disMat:  the user provided dissimilarity matrix, only matters when dis_option
%          is set to 'user'. You can specify disMat as either a full p-by-p matrix,
%          or in upper triangle form. A full dissimilarity matrix must be
%          real and symmetric, and have zeros along the diagonal and non-negative
%          elements everywhere else.  A dissimilarity matrix in upper triangle
%          form must have real, non-negative entries. NaNs are treated as
%          missing values and ignored.  Inf is not accepted.
% isNewSub: i) 1*n vector of 0s or 1s, where
%              1 : the data for the corresponding subject, i.e.,
%                  X(isNewSub == 1, :) are used for prediction only.
% 
%              0 : the data for the corresponding subject, i.e.,
%                  X(isNewSub == 0, :)  are used for estimation only.
% 
%           This option is convenient for computing leave-one-out
%           prediction if desired.
% 
%           ii) If it is a positive integer, say 5, then the last 5 subjects
%              are used for prediction. This option is convenient for obtaining predictions 
%              for a set of new subjects.
% 
%           iii) set isNewSub = [] for the default value, which is no prediction.
% 
% 
% isFPCA : a scalar (0 or 1) that indicates whether functional principal compoenent
%          analysis should be performed on the stringed data. 
%          0 : not perform FPCA
%          1 : perform FPCA and the results are stored in xx. (default)
% isPlot:  a scalar(0 or 1) that indicates whether the stringing function. 
%          0: not plot (default)
%          1: plot
% Details: i) Any unspecified or optional arguments can be set to "[]" for
%             default values;
% =======
% Output:
% =======
% 
%  newXmat: a matrix of n*p contains the stringed data.
%  stringed: the estimated stringing function. newXmat = Xmat(:, stringed);
%  disMatrix: the dissimilarity matrix used in stringing. If not provided
%             by the user, it is calculated based on Xmat(IsNewSub == 1,:); 
%  x      : 1*n cell array of the corresponding random trajectory for each subject identified
%           by the stringing method, where x{i} is the row vector of
%           measurements for the ith subject, i=1,...,n. Those x(IsNewSub
%           == 1) are for prediction. 
%  t_x    : 1*n cell array, t_x{i} is the row vector of time points for the ith
%           subject at which corresponding measurements x{i} are taken,
%           i= 1,...,n. 
%  SubId:  a 0-1 vector indicates which subjects are used for estimation (0) and
%          which are for prediction (1), computed from isNewSub.
%  xx     : the output from FPCA.m based on the subject with isNewSub == 0.
%  pred_x:  1*n cell array, pred_x{i} is the row vector of estimtaed (or predicted)
%            curve for subject i. 
%  pred_pc: n*K matrix of estimated (or predicted) FPC scores, where K is the number of components selected. 
%
%  Note: If isFPCA == 0, the output xx, pred_x and pred_pc are all [].
%         The output x and t_x can be used for further functional data
%         analysis as needed,
%         such as functional principal component analysis and functional regression.
%         See FPCreg, FPCfam and so on.
% 
%  See exampleStringing.m

function [newXmat, stringed, disMatrix, x, t_x, SubId, xx, pred_x, pred_pc] = Stringing(Xmat, dis_option, disMat, isNewSub, isFPCA, isPlot)
   
   [SubId, invalid] = SubForString(Xmat, isNewSub);
   if invalid==1
       return
   end
  
   if isempty(dis_option)
       dis_option = 'euclidean';
   end
   
   if isempty(isPlot)
       isPlot = 0;
   end
   
   if strcmp(dis_option, 'user') == 0     
       disMat = pdist(Xmat(SubId==0, :)', dis_option);
   end
   disMatrix = disMat;   
   n = size(Xmat, 1);
   p = size(Xmat, 2);
   p1=mdscale(disMatrix, 1, 'criterion','metricstress');
   [pos_sort stringed] = sort((p1-min(p1))/range(p1));
   if isPlot ==1 
       a = 1:length(stringed);
       [b ind] = sort(a(stringed));
       figure()
       plot(a,ind,'*-')
       xlabel('Original Order');
       ylabel('UDS Induced Order')
       title('Stringing Function', 'FontWeight','Bold');
   end
   newXmat=Xmat(:,stringed);
   newXtime = repmat(1:p,n,1);
   x = mat2cell(newXmat,repmat(1,1,n),p)';
   t_x =mat2cell(newXtime,repmat(1,1,n),p)';
   
   if isFPCA == 0 
       xx = [];
       pred_x = [];
       pred_pc = [];
   else 
       paramX=setOptions('regular',2, 'bwmu', 3, 'bwxcov', [3, 3], 'selection_k', 'FVE', 'FVE_threshold', 0.99);
       xx = FPCA(x(SubId==0), t_x(SubId==0), paramX);
       [pred_x, pred_pc, xivar] = FPCApred(xx, x, t_x, 2); 
       clear xivar;
   end
   