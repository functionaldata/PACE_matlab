% Example file for how to use function repfFPCA.
%
% The goal is to estimate eigenfunctions and functional principal 
% component scores.

%%% Mortality data, regular dense example

load DataforM2.txt
out1 = 1960:2006;
out2 = 60:100;
m = length(out1);
tn = length(out2);
countries = 1:32;
n = length(countries);
mat = reshape(DataforM2(:,4), m*tn, n)';
 countries = [1,2,4,6,7,8,9,10,11,12,13,14,15,16,17,19,20,21,22,23,24,26,27,28,29,30,32];
 mat = mat(countries, :);
 n = length(countries);
Xmat = reshape(mat, n, tn, m);


verbose = 'off';
param_xi = setOptions('bwmu', 3, 'bwxcov', [3,3], 'selection_k', 'FVE', 'FVE_threshold', 0.95, 'xmu', zeros(1,m), 'verbose', verbose);
ngrid = tn;
design =0;
bwphi = [4,4];
K = 2;
FVE_threshold = 0.9;
mu = [];
xcov = [];
y=[];
t=[];
bwmu = [];
bwxcov = [];
[repf_res] = repfFPCA(design, Xmat, y, t, out1, out2, ngrid, mu, xcov, K, FVE_threshold, param_xi, bwphi, bwmu, bwxcov);
        

%%% sparsified mortality data
load DataforM2.txt
out1 = 1960:2006;
out2 = 60:100;
m = length(out1);
tn = length(out2);
countries = 1:32;
n = length(countries);
mat = reshape(DataforM2(:,4), m*tn, n)';
 countries = [1,2,4,6,7,8,9,10,11,12,13,14,15,16,17,19,20,21,22,23,24,26,27,28,29,30,32];
 mat = mat(countries, :);
 n = length(countries);
Xmat = reshape(mat, n, tn, m);

tmat1 = reshape(repmat(out1, tn, 1), 1, tn*m);
tmat2 = repmat(out2, 1, m);
y= cell(1,n);
t = cell(1,n);
tt = [tmat1; tmat2];
for i = 1:n
    [something,ind] = sort(rand(m,1));
    j = out1(ind(1:32));
    ids = (ismember(tmat1, j)==1);
    t{i} = tt(:, ids);
    y{i} = Xmat(i,ids);
end

verbose = 'off';
param_xi = setOptions('bwmu', 3, 'bwxcov', [3,3], 'selection_k', 'FVE', 'FVE_threshold', 0.95, 'xmu', zeros(1,m), 'verbose', verbose);
ngrid = 20;
design =1;
bwphi = [];
K = 2;
FVE_threshold = 0.9;
mu = [];
xcov = [];
bwmu = [3,3];
bwxcov = [3,3,3];
[repf_res] = repfFPCA(design, Xmat, y, t, out1, out2, ngrid, mu, xcov, K, FVE_threshold, param_xi, bwphi, bwmu, bwxcov);


%%% plot figures for repf_res
phi_all= getVal(repf_res, 'phi_all'); 
xi_pred = getVal(repf_res, 'xi_pred');
psi_all = getVal(repf_res, 'psi_all');

mu =  getVal(repf_res, 'mu');
xi_all = getVal(repf_res, 'xi_all');

varphi = getVal(repf_res, 'varphi');
predysurface = getVal(repf_res, 'predysurface');
out1 = getVal(repf_res, 'out1');
out2 = getVal(repf_res, 'out2');


figure()
subplot(2,2,1)
mesh(out1, out2, phi_all(:,:,1));
title('phi_1')
subplot(2,2,2)
for i = 1:n
    hold all
    plot(out1, xi_pred{1}{i})
end
hold off 
subplot(2,2,3)
mesh(out1,out2, -phi_all(:,:,2));
title('phi_2')
subplot(2,2,4)
for i = 1:n
    hold all
    plot(out1, xi_pred{2}{i})
end
hold off

figure()
subplot(2,2,1)
mesh(out1, out2, mu);
title('$$\mu$$','Interpreter','latex', 'FontSize', 18)

subplot(2,2,2)
k1=1;
k2 = 1;
mesh(out1, out2, varphi(:, :, k2, k1))
title('$$\varphi_{11}$$','Interpreter','latex', 'FontSize', 18)
subplot(2,2,3)
k1 = 1;
k2 = 2;
mesh(out1, out2, varphi(:, :, k2, k1))
title('$$\varphi_{12}$$','Interpreter','latex', 'FontSize', 18)
subplot(2,2,4)
k1 = 2;
k2 = 1;
mesh(out1, out2, -varphi(:, :, k2, k1))
title('$$\varphi_{21}$$','Interpreter','latex', 'FontSize', 18)

