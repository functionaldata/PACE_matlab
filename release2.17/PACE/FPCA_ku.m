%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Estimate FPCA components as in Kneip and Utikal
%
% Input
%
%	t - 1xm grid of input values
%   f - nxm matrix of function values - ith column corresponds to input t(i)
%   w - optional 1xm vector of weights for inner product between functions
%   (default is even weights)
%	selection_k - If 'FVE', then use FVE_threshold as criteria for selecting number of principal components
%	If a positive number : user-sepcified number of principal components
%  FVE_threshold - default 0.95 
% Output
%	refer to the output of PCA

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [no_opt,sigma,lambda,phi,eigen,xi_est,xi_var,mu,muDense,bw_mu,xcov,bw_xcov,xcovfit,AIC,BIC,FVE,y_pred,y_predOrig,y_predDense,out1,out21,y,t, regular, rho_opt, sigmanew, mucopy, phicopy, eigencopy, out1copy, out21copy, xcovcopy, xcovfitcopy, xcorr, ops] = FPCA_ku(t, f, selection_k,FVE_threshold);

	if(nargin < 4)
		FVE_threshold = 0.95;
	end
	
	if(nargin < 3)
		selection_k = 2;
	end
	
	default_FVE = 0; % the switch of FVE method.
	if (isnumeric(selection_k))
		K = selection_k;
	elseif(strcmp(selection_k,'FVE'))
		default_FVE = 1;
	else
		fprintf('When regular is 3, Only FVE approach is available for selecting principal components\n');
		default_FVE = 1;
	end
	
	n = size(f, 1); m = length(t);
	mu = mean(f);
	M = zeros(n);
	
	for i = 1:n
		for j = 1:n
			M(i, j) = trapz(t, (f(i, :) - mu).*(f(j, :) - mu));
		end
	end
	
	[p lam] = eig(M);
	[lam, ord] = sort(diag(lam), 'descend');
	p = p(:, ord);
	cumlam = cumsum(lam)/sum(lam);
	ind_k = find(cumlam>FVE_threshold);
	ind_k = ind_k(1);  % the smallest # of components  
	if  default_FVE == 1
		K = ind_k;  % if user doesn't specify K
	end
	
	xi = p(:, 1:K)*sqrt(diag(lam(1:K)));
	phi = zeros(K, m);
	lam = lam /n; % divided by subject number 
	
	%transform back the eigenfunction corresponding to original data
	for r = 1:K
		phi(r, :) = sum(diag(xi(:, r))*f)/sum(xi(:, r).^2);
	end
	
	%estimate the extra output for PACE
	no_opt = K;
	sigma = [];
	lambda = lam(1:K)';
	phi = phi';
	eigen = [];
	xi_est  = xi;
	xi_var = [];
	mu;
	muDense = [];
	bw_mu = [];
	xcov = [];
	bw_xcov = [];
	xcovfit = [];
	AIC = [];
	BIC = [];
	FVE = cumlam;
	y_pred = [];
	y_predOrig = repmat(mu,n,1)+ xi_est * phi';
	y_predOrig = mat2cell(y_predOrig,ones(1,n),m)';
	y_predDense = [];
	out1 = [];
	out21 =[];
	y = mat2cell(f,ones(1,n),m)';
	t = mat2cell(repmat(t,n,1),ones(1,n),m)';
	regular = 3; 
	rho_opt = [];
	sigmanew = []; 
	mucopy= []; 
	phicopy = [];
	eigencopy = [];
	out1copy= [];
	out21copy= [];
	xcovcopy= [];
	xcovfitcopy = [];
	xcorr = [];
	ops =[];

end
