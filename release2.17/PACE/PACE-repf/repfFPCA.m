%% 
% This is the main function for repeated functional data modeling. 
% Input:
% design: 0 regular dense design
%         1 irregular or sparse design
% Xmat: the data array n*tn*m for design =0, set it to [] for design = 1;
% y: 1*n cell, y{i} values corresponding to t{i} for design =1, set to [] for design =0;
% t: 1*n cell, t{i} is a 2*ni matrix, contains all pairs (s,t) for subject
%      i, set to [] for design =0;
% out1: the grid for s, length(out1) = m, for design =0;  
%       for design =1, out1 = unique(tt(1,:)), tt = cell2mat(t);
% out2: the grid for t, length(out2) = tn for design =0;
%       for design =0, out2 = unique(tt(2,:)), tt= cell2mat(t);
% ngrid: the number of grids for t, when estimating xcov.
% mu: length(out2)*length(out1) matrix, the user provided mean surface.
%     If it is [], then empirical esitmator will be used for design =0, and
%     smoothing estimator for design =1;
% xcov: ngrid*ngrid*length(out1), the user provided covaraince surfaces.
%       If it is [], then empirical estimator will be used for design =0,
%       and smoothing estimator for design =1;
% K: the user provided number of components for the first step FPCA
% FVE_threshold: for choosing the number of components, if K is not provided
% param_xi: the parameters for the second FPCA
% bwphi: the bandwidth for the addtional smoothing for phi(t|s), in the case that empirical cov is used. 
%        If no addtional smoothing is wanted, set it to []. 
% bwmu: the bandwidth for smoothing mu; set to [] if not needed.
% bwxcov: the bandwith for smoothing xcov; set to [] if not needed.

% Output:
% K: the number of FPC components chosen for the first step FPCA
% pk: the number of FPC components chosen for the second step FPCA based on
%     the working process $\xi_k(s)$.
% P:  the max of pk.
% FVEk: the FVE explained. 
% out1: the grid for s.
% out2: the grid for t.
% mu: the estimated mean function of X(t,s)
% xcov: the estimated covariance function, cov(X(t1,s), X(t2,s))
% xi_all(n, m, K): the FPC scores from the first step FPCA, used as working
%                  processes for the second step FPCA. 
% xi_pred: the predicted values of xi_all, after fitting the second step FPCA. 
% phi_all: the eigen functions from the first step FPCA, 
%          phi_all(.,s, k) is the eigen functions 
% lambda_all: the eigen values from the first step FPCA, 
%          lambda_all(s,k) is the kth eigen value from FPCA of $X_i(., s)$.
% zeta(n, P,K): the eigen value from the second step FPCA, 
% psi_all(m, P,K): the eigenfunctions from the second step FPCA
% predysurface: the predicted X_i(t,s).

% Reference: Chen K and Mueller HG (2012) Modeling repeated functional
% observatsion.

function [repf_res] = repfFPCA(design, Xmat, y, t,out1, out2, ngrid, mu, xcov, K, FVE_threshold, param_xi,bwphi, bwmu, bwxcov)
if isempty(FVE_threshold)
    FVE_threshold = 0.85;
end

% if the design is dense regular, compute empirical mean and
% covariance G(t1,t2|s).
if design == 0 
   n = size(Xmat,1);
   tn = size(Xmat, 2);
   m = size(Xmat, 3);
   ngrid = tn;
   if isempty(mu)
       mu = reshape(mean(Xmat), tn, m);
   end
   if isempty(xcov)
       xcov = zeros(ngrid, ngrid, m);
       for i = 1:m
            xcov(:,:,i) = (Xmat(:,:,i)-repmat(mu(:,i)', n, 1))'*(Xmat(:,:,i)-repmat(mu(:,i)', n, 1));
       end
   end
   tmat1 = reshape(repmat(out1, tn, 1), 1, tn*m);
   tmat2 = repmat(out2, 1, m);
   tmat = repmat([tmat1; tmat2], n, 1);
   t = mat2cell(tmat, repmat(2, n, 1), tn*m)';
   mat = reshape(Xmat, n, tn*m);
   y = mat2cell(mat, repmat(1, n, 1), tn*m)';
end

% if the design is random, compute mean and 
% covariance G(t1,t2|s) through 2-d and 3-d smoothing. 
if design == 1
   if isempty(mu) 
        tt = cell2mat(t);  
        yy = cell2mat(y);  
        fprintf(1,'Obtain smoothed mean surface\n');
        rawmean = getRaw2dMean(tt, yy);
        outmu = GetMesh(out1, out2);
        win1 = ones(1,size(rawmean.tt,2));
        [invalid, mu] = mean2d(bwmu,'epan',rawmean.tt,rawmean.yy',win1,outmu,rawmean.count);
        if invalid == 1 
           return
        end
        mu = reshape(mu, length(out2), length(out1));
        clear rawmean win1 tt yy; 
   end
   if isempty(xcov)
       error = 1;
       fprintf(1,'Obtain smoothed covariance function\n');
       out21 = linspace(min(out2), max(out2), ngrid);
       rcov = getRawRepfCov(y, t, mu, out1, out2);           %obtain raw covariance;
       if error == 1
         tpairn = rcov.tpairn;
         tneq = find(tpairn(1,:)~=tpairn(2,:));
         rcov.cyy = rcov.cyy(tneq);
         rcov.tpairn = rcov.tpairn(:,tneq);
         rcov.count = rcov.count(tneq);
       end
       win=ones(1,length(rcov.cyy));
       [invalid,xcov]=GetRepfCov(bwxcov,'epan',rcov.tpairn,rcov.cyy',win,out21, out1,rcov.count);  %smooth raw covariance;
       clear tpairn win rcov out21;
       if invalid ==1 
         return
       end
   end
   n = length(y); 
   m = length(out1);
   tn = length(out2);
end

xx = cell(1, m);
kxx = zeros(m, 1);
FVEk = zeros(m, 1);

% The first step FPCA

for i = 1: m 
    xx{i} = FPCAstep1(out2,K, ngrid, xcov(:,:,i));
    kxx(i) = getVal(xx{i}, 'no_opt');
    FVEvalue = getVal(xx{i}, 'FVE');
    FVEk(i) = find(FVEvalue >= FVE_threshold, 1);    
end
if isempty(K)
    K = max(FVEk);
elseif (K > min(kxx))
    K = max(FVEk);
    fprintf('User defined K is too big. Reset to K that minimum FVE = 0.85 over all s \n')
end


% Determine the order of eigenfunctions according to smoothness. 
phi_all = zeros(tn, m, K);
lambda_all = zeros(m, K);
for i = 1:m
    phi = getVal(xx{i}, 'phi');
    lambda = getVal(xx{i}, 'lambda'); 
    if i == 1
        for k = 1:K
            phi_all(:, i, k)= phi(:,k);
            lambda_all(i,k) = lambda(k);
        end
    else
        set = 1:K;
        for k = 1: K
            lk = set(1);
            lksign = 1;
            for l = set 
                bestsign = 1;
                if (norm(phi(:,l) - phi_all(:, i-1, k)) > norm(phi(:,l) + phi_all(:, i-1, k)))
                    bestsign = -1;
                    if (l == lk) 
                        lksign = -1;
                    end
                end
                if (norm(bestsign*phi(:,l) - phi_all(:, i-1, k)) < norm(lksign*phi(:,lk) - phi_all(:, i-1, k)))
                    lk = l;
                    lksign = bestsign;
                end
            end
              phi_all(:, i, k)= lksign * phi(:, lk); 
              lambda_all(i,k) = lambda(lk);
              set = setdiff(set, lk);
        end
    end
end

% Additional smoothing on eigenfunctions if needed.
if ~isempty(bwphi)
  kernel = 'epan';
  xin = GetMesh(out1, out2);
  win =  ones(1, m*tn);
  for k = 1:K
       yin = reshape(phi_all(:,:,k), tn*m, 1);
       [invalid phi_all(:,:,k)] = mullwlsk_M(bwphi, kernel, xin, yin, win, out1, out2);
  end
end

% Compute the FPF scores, which will be used as working process for the
% second step FPCA
xi_all = NaN(n, m, K);
for j= 1:m
   [x, t_x insub] = GetXT(j, n, out1, y, t);
   ii = 0;
   for i = 1:n
       if insub(i) == 1
           ii = ii+1;
           xinterp  = interp1(t_x{ii}, x{ii}, out2 ,'spline');
           for k = 1:K             
               prod = (xinterp'-mu(:,j)).* phi_all(:, j, k);
               xi_all(i, j, k) = trapz(out2, prod);
           end
       end
   end
end

x_xi = cell(K, n);
t_xi = cell(K, n);
for k  = 1:K
  for i = 1:n
   idm = ~isnan(xi_all(i,:,k));
   x_xi{k, i} = xi_all(i,idm, k);
   t_xi{k, i} = out1(idm);
  end
end
 
% Second step FPCA
xixi = cell(1, K);
pk = zeros(K,1);
for k = 1: K
  xixi{k} = FPCA(x_xi(k,:), t_xi(k,:), param_xi);
  pk(k) = getVal(xixi{k}, 'no_opt');
end
P = max(pk);
zeta = zeros(n, P, K);
psi_all= zeros(m, P, K);
xi_pred = cell(1, K);
for k= 1:K 
  psi_all(:, 1:pk(k), k)= getVal(xixi{k}, 'phi');
  zeta(:, 1:pk(k), k) = getVal(xixi{k}, 'xi_est');
  xi_pred{k} = getVal(xixi{k}, 'y_pred');
end

varphi = zeros(tn, m, P, K);
predysurface = mat2cell(repmat(mu, n, 1), repmat(tn, 1, n), m)' ;
%adjust = repmat(mu2(:,1)', tn, 1).*phi_all(:,:,1) + repmat(mu2(:,2)', tn, 1).*phi_all(:,:,2);
adjust = 0;

for k1 = 1:K
    for k2 = 1:P
        varphi(:, :, k2, k1) = phi_all(:, :, k1) .* repmat(psi_all(:, k2, k1)', tn, 1);
        for i = 1:n
            predysurface{i} = predysurface{i} + zeta(i, k2,k1)* varphi(:,:,k2,k1) + adjust;
        end
    end
end

Xnames = {'K', 'P', 'pk','FVEk', 'out1', 'out2', 'mu', 'xcov', 'xi_all', 'xi_pred', 'phi_all', 'lambda_all', 'zeta', 'psi_all', 'varphi', 'predysurface',  'names'};
repf_res = {K, P, pk, FVEk, out1, out2, mu, xcov, xi_all, xi_pred, phi_all, lambda_all, zeta, psi_all, varphi, predysurface, Xnames};

end