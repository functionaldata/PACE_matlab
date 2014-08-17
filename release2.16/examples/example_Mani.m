addpath(genpath('~/fda/release2.14/'));
addpath(genpath('~/fda/codes/'));



 
% generate data: regularly spaced, no observation error
N = 200;
X = cell(1,N);
T = cell(1,N);
R2 = ones(1,N);

alph = normrnd(0,.75,N,1) ;
cT = (-125:5:125)/50;
eps = normrnd(0,.5,1,N);
for i=1:N
    X{i} = normpdf(cT,alph(i),.1);
    T{i} = cT;
    R2(i) = alph(i)+ eps(i);
end
Y = R2;
t = cT;




p = setOptions( 'regular', 2, 'selection_k', 'FVE','FVE_threshold', 0.9,'screePlot',1, 'designPlot',1,'corrPlot',1,'numBins',0, 'verbose','on');

[Y,Outliers,S] = maniMDS(X,T,3,2,[],[],[],[],[],[],p);


