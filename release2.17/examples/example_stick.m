% Example file for how to use stick.m.
% Sparse, irregular, random design, 80 subjects, 21 measurements on each subject.
% Curves are generated each by 3 FPCs generated from independent Gaussian, with
%     s.d 3,2 and 1, respectively.
% Gaussian noise is added at each measurement with s.d 0.5.

% The goal is to calculate the functional correlation coefficient (named 'stickiness 
%    coefficient' )of two sparse, irregular observed realizations of functional data.

% Reference: 
%     Gottlieb, A. and M\"{u}ller, H.G. (2012). A Stickiness Coefficient 
%     for Longitudinal Data. Computational Statistics and Data Analysis (Accepted).


p = path;
isExist = regexp(p, 'PACE');
if isempty(isExist) == 1
addpath(genpath('../PACE/'));
end


% time domain is [0,10] 
% all the observations are taken equidistantly.


%generate data set
%clear all;


rng(123, 'twister');
ncohort = 80;
nobs = 50;
lint_x = 10; 
numeig_x = 3; 

tx = 0:lint_x/nobs:lint_x;           % observation time points 

phi_x = xeig(tx,lint_x,numeig_x);    % 3 eigenfunctions

lambda_x = [3 2 1].^2;
pc_x = randn(ncohort,3)*diag(sqrt(lambda_x));         % 3 FPC score 
mu_x = mu_true(tx,lint_x);          % true mean function 




x = cell(1,ncohort); 
t = mat2cell(repmat(tx,ncohort,1),ones(1,ncohort))';


for i=1:ncohort
    %generate the repeated measurements  with measurement errors N(0,0.25)
    x{i} = mu_x+pc_x(i,:)*phi_x+0.5*randn(1,length(tx));
end


% specify parameters for implementing FPCA analysis on the functional data.
param_X = setOptions('selection_k', 'FVE','FVE_threshold', 0.9);
% param_X=[];

% Implement stickiness estimation
[stickiness_coefficient] = stick(x, t, param_X)
%
