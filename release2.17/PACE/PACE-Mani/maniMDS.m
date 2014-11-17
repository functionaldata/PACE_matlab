%=========== Description: ===========
%
% Main program to perform nonlinear dimension reduction via ISOMAP or
% Penalized-ISOMAP, which use geodesic distances and multidimensionl
% scaling techniques (maniMDS means manifold multidimensional scaling).
% 
% Reference: Chen, D., Müller, H.G. (2012). Nonlinear manifold
% representations for functional data. Annals of Statistics (in press)
%
%================ Input Arguments: ================
%
%Input X: 	1*N cell array; X{i} is a 1*D_i vector with measurements for the ith
% 		subject, i=1,...,N. 
%
%Input	T: 	1*N cell array; T{i} is a 1*D_i vector with the time location
% 		measurements for the ith subject, i=1,...,N. 
%
%Input	datatype: integer 0-4, determing how the data is treated (dense/sparse)
% 		0 - sparse, use PACE and special distance (Peng and M\"uller 2008)
% 		1 - sparse, use PACE and L2 distance [default] 
%		2 - regular/with missing data, use PACE and L2 distance
%		3 - dense with noise/missing values, use local linear smoothing
%		4 - dense and balanced, without noise
%       5 - datatype 4, but use custom distance
%
%Input d: 	integer, output dimension [default 3]
%
%Input t: 	vector, output time grid [default balanced with 30 points]
%
%Input method: string giving the method name
% 		'isomap' - use ISOMAP (Tenenbaum, 2000)
% 		'p-isomap' - penalized isomap
%
%Input kernel: string giving the kernel name
% 		'quar' - quartic kernel
% 		'epan' - epanechikov kernel
% 		'gauss' - gaussian kernel
% 		'gausvar' - variant of gaussian kernel
% 		'rect' - rectangular kernel
%
%Input K:	integer indicating how many neighbors are used to find geodesic
% 		distance
%
%Input h:	kernel bandwidth used to predict functional trajectories based on FMC
% 
%Input delta: truncation parameter used in penalized isomap.
%
%Input pace:	parameter used in PACE, see FPCA
%
%Input Dis: If input datatype == 5, a custom distance matrix
%
%================ Output Arguments: ================
%
%Output Y:		estimated functional manifold components
% 
%Output Outliers:	0/1 indicator; 1 if the corresponding subject is outlier
%
%Output S:	 	a structure that contains Y, Outliers, Dis, t, N, M, X, T, X_pace,
%			X_reg, datatype, d, method, K, h, delta, mse, kernel, FVE, mu, xi,
%			phi, lambda.
% 
%			1) Y - N*d matrix of estimated functional manifold components
%
%			2) Outliers - 1*N logical indicating whether (1) or not (0) each 
%					subject is an outlier
%
%			3) Dis - N*N matrix of estimated L2 distances between trajectories
%
%			4) t - 1*M row vector containing the output time grid
%            
%			5) N - integer, number of subjects
%
%			6) M - integer, number of timepoints in output time grid t
%
%			7) X - see the input of the same name
%
%			8) T - see the input of the same name
%
%			9) X_pace - cell obtained from FPCA(X,T,pace);
%
%			10) X_reg - N*M matrix containing estimated data trajectories 
%					evaluated at regular time points
%
%			11) datatype - see the input of the same name
%
%			12) d - see the input of the same name
%
%			13) method - see the input of the same name
%
%			14) K - see the input of the same name
%
%			15) h - see the input of the same name
%
%			16) delta - see the input of the same name
%
%			17) mse - array containing estimated mean squared errors,
%				 	obtained by 10-fold CV, for candidate values of 
%					K, h, and delta. If these inputs each are set to [],
%					the default results in 5*5*10 array (see par10cv.m).
%
%			18) kernel - see the input of the same name
%
%			19) FVE - numeric fraction of variance explained (more accurately, 
%					fraction of dimension explained)
%
%			20) mu - 1*M vector, mean trajectory corresponding to output X_reg
%
%			21) xi - N*d matrix, predictions for first d random PC coefficients, 
%					obtained from output X_reg
%
%			22) phi - M*d matrix, estimated eigenfunctions, obtained from output X_reg,
%					corresponding to the largest d eigenvalues
%
%			23) lambda - 1*d matrix, first d estimated eigenvalues (variances of PC scores)
%					obtained from output X_reg
%
%%%% see also PCA.m
%







function [Y,Outliers,S] = maniMDS(X,T,datatype,d,t,method,kernel,K,h,delta,pace,Dis)

if nargin<12 Dis = []; end
if nargin<11 pace = []; end
if nargin<10 delta = []; end
if nargin<9 h = []; end
if nargin<8 K = []; end
if nargin<7|isempty(kernel) kernel = 'epan'; end
if nargin<6|isempty(method) method = 'p-isomap'; end
if nargin<5 t = []; end
if nargin<4|isempty(d) d = 3; end
if nargin<3|isempty(datatype) datatype = 1; end
if isempty(t)&~isempty(getVal(pace,'newdata')) t = getVal(pace,'newdata'); end

N = length(X);
t_pooled = unique(cell2mat(T));
if datatype==5|datatype==4|datatype==3
    if isempty(t) t = t_pooled; end
    M = length(t);
    if datatype==5|datatype==4
        X_reg = reshape(cell2mat(X),M,N)';
    else
        X_reg = lwlsCV(X,T,t);
    end
    % create X_cv on grid t_pooled for Cross-validation
    if length(t)~=length(t_pooled)|~all(t==t_pooled)
        X_cv = lwlsCV(X,T,t_pooled);
    else
        X_cv = X_reg;
    end
    if datatype==4|datatype==3
        Dis = L2_distance(X_reg',X_reg',1)*sqrt(range(t)/(M-1));
    end
    X_pace={};
else
    if isempty(t) t = min(t_pooled):range(t_pooled)/29:max(t_pooled); end
    M = length(t);
    if isempty(pace)
        if datatype==2
            regular = 2;
        else
            regular = 0;
        end
        pace = setOptions('bwmu',0.45517,'bwxcov',[0.18207,0.18207],'regular',regular,'selection_k','FVE','FVE_threshold',0.99,'newdata',t,'verbose','off');
    end
    X_pace = FPCA(X,T, setfield(pace,'newdata',t) );
    if datatype == 0
        Dis = spadis(X_pace);
    else
        xi_est = getVal(X_pace,'xi_est');
        Dis = L2_distance(xi_est',xi_est',1);
    end
    X_reg = reshape(cell2mat(getVal(X_pace,'y_pred')'),M,N)';
    % create X_cv on grid t_pooled for Cross-validation
    if (length(t)~=length(t_pooled)|~all(t==t_pooled))&(isempty(K)|isempty(delta)|isempty(h))
        regular = getVal(pace,'regular');
        tmptpace = setOptions('regular',regular,'selection_k','FVE','FVE_threshold',0.95,'newdata',t_pooled,'verbose','off');
        tmptX_pace = FPCA(X,T,tmptpace);
        X_cv = reshape(cell2mat(getVal(tmptX_pace,'y_pred')'),length(t_pooled),N)';
        clear tmptpace tmptX_pace;
    else
        X_cv = X_reg;
    end
end
Dis = (Dis+Dis')/2;
if strcmp(method,'isomap') delta = 0; end
if isempty(K)|isempty(delta)|isempty(h)
    [Kcv,deltacv,hcv,mse] = par10cv(X,T,datatype,X_cv,t_pooled,kernel,Dis,method,d,K,delta,h);
else 
    Kcv = K;
    hcv = h;
    deltacv = delta;
    mse = NaN;
end
Outliers = logical(zeros(1,N));
[Y,Index,FVE] = maniMethods(Dis,method,d,Kcv,deltacv);
if length(Index)<N
    Dis = Dis(Index,Index);
    Outliers(Index) = logical(ones(1,length(Index)));
    Outliers = ~Outliers;
    N = size(Dis,1);
end

if nargout==3
    S.Y = Y;
    S.Outliers = Outliers;
    S.Dis = Dis;
    S.t = t;
    S.N = N;
    S.M = M;
    S.X = X;
    S.T = T;
    S.X_pace = X_pace;
    S.X_reg = X_reg(~Outliers,:);
    S.datatype = datatype;
    S.d = d;
    S.method = method;
    S.K = Kcv;
    S.h = hcv;
    S.delta = deltacv;
    S.mse = mse;
    S.kernel = kernel;
    S.FVE = FVE;
    if datatype>=3
        S.mu = mean(S.X_reg);
        [xi,phi,lambda] = pceig(X_reg,t,d);
        S.xi = xi;
        S.phi = phi';
        S.lambda = lambda;
    else
        S.mu = getVal(X_pace,'mu');
        tmpt_xi = getVal(X_pace,'xi_est');
        S.xi = tmpt_xi(~Outliers,:);
        S.phi = getVal(X_pace,'phi');
        S.lambda = getVal(X_pace,'lambda');
    end
end

end

