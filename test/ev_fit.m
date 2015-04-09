% 
tic
% no_eig = 5;
K_true = 50;
M = 1;
n_train = 50;
n_test = 0;
n_total = n_train + n_test;
shift = 0;
lambda_true = (1:K_true).^(-2);
seed = 1;
PC_dist = 'norm';
sigma2_true = 0.1; % ASSUME observation error
K_max = 10;
unif_on = 4:2:10;

grid_pts = 0:0.02:1;
[basis_true, ~] = makeBasis(grid_pts, K_true, 0);

if (ispc)
    addpath 'D:\Documents\My Documents\Dropbox\Research\PACE_matlab\release2.16\PACE'
    addpath 'D:\Documents\My Documents\Dropbox\Research\PACE_matlab\release2.16\PACE\PACE-FAM'
end
if (isunix)
    addpath '~/PACE_matlab/release2.17/PACE'
    addpath '~/PACE_matlab/release2.17/PACE/PACE-FAM'
end

rng(seed);
for (mc = 1:M)
    mc
    [Y, group] = sim_obs(n_total, 1/2, lambda_true, lambda_true, ...
        basis_true, basis_true, PC_dist, sigma2_true);
   
    dat_all = num2cell(Y, 2)';
    t_all = num2cell(ones(n_total ,1) * grid_pts, 2)';
    
    [dat_thin, t_thin] = thinObs(dat_all, t_all, unif_on);
    
    regular = 0; kernel = 'epan'; error_type = 1; rho = 0; verbose = 'off'; 
    method = 'CE'; shrink = 0; 
    out1_t = grid_pts;
    p1=setOptions('selection_k', K_max, 'regular', regular, 'method', method, 'kernel', ...
    kernel, 'numBins', 0, 'newdata', out1_t, 'error', error_type, 'screePlot', 0, ...
    'rho', rho, 'verbose', verbose, 'ls_fit', 1);
    res1 = FPCA(dat_all, t_all, p1);
    
    p0=setOptions('regular', regular, 'method', method, 'kernel', ...
    kernel, 'numBins', 0, 'newdata', out1_t, 'error', error_type, 'screePlot', 0, ...
    'rho', rho, 'verbose', verbose, 'ls_fit', 1, 'newdata', [0.02, 0.98]);
    res0 = FPCA(dat_all, t_all, p0); 
    
end
time_elapsed = toc;

