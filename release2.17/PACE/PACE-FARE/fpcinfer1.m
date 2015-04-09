function [cbmu,cbfhat,cbcov,cbphi,btneig] = fpcinfer1(x,res,alpha,B,kern,xrange,...
    nbins,nbins1,ngrid,ngrid1)
%   alpha:  level of confidence bands, optional; 100(1-alpha)% confidence
%           bands will be created.  alpha is .05 by default.     

id = getVal(res,'id');
x = x(id);
n = length(x);              
out1 = getVal(res,'out1');
out21 = getVal(res,'out21');
bw = getVal(res,'bw');
noeig = getVal(res,'noeig');

if isempty(alpha)
    alpha = 0.05;
end
if isempty(B)
    B = 200;
end
if isempty(kern)
    kern = 'gauss';
end
if isempty(xrange)
    xrange = [min(out1),max(out1)];
end
if isempty(nbins)
    nbins = 101;
end
if isempty(nbins1)
    nbins1 = 101;
end
if isempty(ngrid)
    ngrid = 51;
end
if isempty(ngrid1)
    ngrid1 = 30;
end

for i = 1:n
    x{i} = x{i}(x{i}<=xrange(2) & x{i}>=xrange(1));
end

btmu = zeros(B,length(out1));
btcov = zeros(B,length(out21),length(out21));
btphi = zeros(B,length(out1),noeig);
btfhat = zeros(B,length(out1),n);
btneig = zeros(B,1);
for i = 1:B
    btx = mysample(x,n,1);
    
    % Estimate the mean density function
    btxx = cell2mat(btx);
    btmu(i,:) = hist2fx(btxx,nbins,[],bw.bw1,kern,out1,xrange); 
    
    % Estimate the covariance surface and eigenvectors
    x2Dtmp = cell(1,n);
    for j = 1:n
        x2Dtmp{j} = nchoosek(btx{j},2)';
    end
    x2D = cell2mat(x2Dtmp);
    
    if bw.covmethod == 0
        f2D = hist2fx2D(btx,x2D,ngrid,[],bw.bw2,kern,ngrid1,out21,xrange);
        fmu21 = interp1(out1,btmu(i,:),out21,'spline');
        Ghat = f2D-fmu21'*fmu21;
    else
        [mutmp,outtmp] = hist1D(btxx,nbins,[],xrange);
        f2dtmp = hist2D(x2D,nbins,[],xrange);
        Ghattmp = f2dtmp-mutmp'*mutmp;
        yin3 = Ghattmp(:);
        [xin31,xin32] = meshgrid(outtmp);
        xin3 = [xin32(:),xin31(:)]';
        clear xin31 xin32;
        indneq = find(xin3(1,:)~=xin3(2,:));
        yin3 = yin3(indneq);
        xin3 = xin3(:,indneq);
        win3 = ones(1,length(yin3));
        [invalid,Ghat] = mullwlsk(bw.bw2,kern,xin3,yin3,win3,out21,out21);
    end
    Ghat = (Ghat+Ghat')/2;
    
    [lambda,phitmp,btcov(i,:,:)] = geteigen(Ghat,out1,out21,1);
    noeigtmp = size(phitmp,2);
    btphi(i,:,:) = phitmp(:,1:noeig);
    
    % Estimate the density function f for each subject
    btx1 = cell(1,n);
    for j = 1:n
        btx1{j} = sort(mysample(x{j},length(x{j}),1));
    end
    [btfhat(i,:,:),xihat,btneig(i)] = getfhat(btx1,xrange,btmu(i,:),phitmp,out1,0,...
        noeigtmp,kern,nbins1);
end

cbmu = quantile(btmu,[alpha/2,1-alpha/2],1);
cbcov = quantile(btcov,[alpha/2,1-alpha/2],1);
cbphi = quantile(btphi,[alpha/2,1-alpha/2],1);
cbfhat = quantile(btfhat,[alpha/2,1-alpha/2],1);

% mu = getVal(res,'fmuhat');
% covfit = getVal(res,'covfit');
% phi = getVal(res,'phi');
% fhat = getVal(res,'fhat');

