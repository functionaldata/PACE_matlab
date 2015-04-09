function [nc_kcfc,initCluster,reCluster,newCluster,idselect,SSE,groupM,groupSSE] = kCFC(y,t,clustopt)
% KCFC
%
% Last modified: 2015/3/5 by Pai-Ling Li
%
    %  1. Set parameters of clustering method for initial clustering
    nc_in = clustopt.nc;
    P = clustopt.P;
    iseed = clustopt.seed;
    initD = clustopt.initD;

    % for reclassification
    M = clustopt.M;
    iter = clustopt.iter;
    moveprop = clustopt.moveprop;
    op_cv = clustopt.op_cv;

    %  2. FPCA  (all curves)      
    %yy = FPCA(y,t,P);
    yy = FPCAw(y,t,P);
    Fpcs = getVal(yy, 'xi_est');
    P_out = getVal(yy, 'no_opt');
    FVE = getVal(yy, 'FVE');
    display(['fclust_kcfc - dimension of FPCs for initial clustering = ', num2str(P_out)]);
    disp(['prop of eigenvalues = ',num2str(FVE(P_out))]);

    %  3. Initial clustering: clustering the FPCs by K-means algorithm
    rand('state',iseed);
    [idc,C,sumd,D] = kmeans(Fpcs(:,[1:P_out]),nc_in,'distance',initD);
    initCluster = idc;

    %  4. Iterative clustering : 
    [reCluster,newCluster,idselect,groupSSE,groupM] = iter_kCFC(y,t,M,initCluster,iter,op_cv,moveprop);

    SSE = groupSSE(idselect,:);
    M_out = groupM(idselect,:);

    [unique_c] = unique(newCluster);
    nc_kcfc = length(unique_c);
    if nc_kcfc < nc_in
       disp(['KCFC - ',num2str(nc_in),' clusters are combined to ', num2str(nc_kcfc),' clusters.']);
    end
end
