% Example: Berkeley Growth Data
% Method: K-center Functional Data Clustering (KCFC) 
%
% Last modified: 2015/3/5 by Pai-Ling Li, before version controlled by git

% clear all;%clc
% addpath(genpath('69p4Chiou'))
% addpath(genpath('kcfc-PACE215'))
addpath(genpath('../PACE/'));

% -------------------------------------------------------------------------
% I n p u t  D a t a 
% -------------------------------------------------------------------------
nc = 2;                 
file_in = 'growth.dat';          
file_out = 'growth.out';     
nrow = 2883;  
[idsubj,idclass,data] = datap(file_in,nrow); 

% -------------------------------------------------------------------------
% D a t a   C l u s t e r i n g   a n a l y s i s
% -------------------------------------------------------------------------

% Set options for inital clustering
clustopt.nc = nc;            % number of cluster
%method = 'CE';
method = 'IN';               % method for estimating FPC scores
clustopt.P = setopts('selection_k', 'FVE', 'FVE_threshold', 0.90, 'regular', 2, 'numBins', 0, 'verbose', 'off', 'method', method, 'ngrid', 288, 'weight', 1, 'error', 0);
clustopt.seed = 232316;            % seed for kmeans algorithm.
clustopt.initD = 'sqEuclidean';    % distance measure of k-means: 'sqEuclidean', 'cityblock', 'cosine', 'correlation', 'Hamming'

% Set optinos for reclustering
clustopt.M = setopts('selection_k', 'FVE', 'FVE_threshold', 0.80, 'regular', 2, 'numBins', 0, 'verbose', 'off', 'method', method, 'ngrid', 288, 'weight', 1, 'error', 0);
clustopt.iter = 50;         % maximum number of iteration
clustopt.moveprop = 0;      % real number between 0 and 1. 
                            % stop reclassification if the number of moving curves less than moveprop*n 
clustopt.op_cv = 0;         % = 0, perform kcfc without leave-one-curve-out cross-validation
                            % = 1, perform kcfc with leave-one-curve-out cross-validation (time comusing)
                            
% Clustering by KCFC
y = mat2cell(data.Yin,ones(93,1))';
t = mat2cell(data.Tin,ones(93,1))';
disp(['nc = ',num2str(clustopt.nc)])
[nc_kcfc,initCluster,reCluster,newCluster,idselect,SSE,groupM,groupSSE] = kCFC(y,t,clustopt);

% Comparing the clustering results with gender classification: 
% Adjusted Rand Index (aRand) and Correct Classification Rate (cRate)
[aRandfpca_init,cRatefpca_init,orderfpca_init,tablefpca_init] = cmp2p(initCluster,idclass) % for initial clustering
[aRandfpca_new,cRatefpca_new,orderfpca_new,tablefpca_new] = cmp2p(newCluster,idclass)      % for reclustering


% Plots
idx1 = find(newCluster == 1);
idx2 = find(newCluster == 2);

figure; hold on
for i = 1:length(idx1)
    if idclass(idx1(i)) == 1
        plot(t{idx1(i)},y{idx1(i)},'b')
    elseif idclass(idx1(i)) == 2
        plot(t{idx1(i)},y{idx1(i)},'r')
    end
end
set(gca,'fontsize',20,'ylim',[60 200],'xlim',[1 18])
%print('latex\growth_fig\cluster1','-depsc')
%close all

figure; hold on
for i = 1:length(idx2)
    if idclass(idx2(i)) == 1
        plot(t{idx2(i)},y{idx2(i)},'b')
    elseif idclass(idx2(i)) == 2
        plot(t{idx2(i)},y{idx2(i)},'r')
    end
end
set(gca,'fontsize',20,'ylim',[60 200],'xlim',[1 18])
%print('latex\growth_fig\cluster2','-depsc')
%close all

