% Example for Stringing.m 

% The p dimensional predictor Xmat is generated from gaussian distribution with 
% a random covariance matrix KK. 
% Response Y is generated from Y = X\beta + error, where beta is generated from
% a gaussian density and is scrambled to random order.   
%
% 60 subjects for estimation; 50 subjects for prediction. 
% The dimension p = 100.

% We perform stringing on Xmat first to get stringed x, and then use functional 
% regression to predict Y as \int x(t)\beta(t).  The mean squared error are computed. 

ntrain =60;
n= ntrain+50;
p=100;

x = 1: p;
beta = exp(-(x-p/2).^2/(2*30^2));
beta = beta';
new = mysample(1:p, p , 0);
beta = beta(new);

Kraw = rand(p);
Kraw = (Kraw + Kraw')/2;
for i = 1 : p 
    Kraw(i,i) = 1; 
end
[eigen d] = eigs(Kraw, p-1,'lm');                                  
d(d <= 0)= 0; 
KK = eigen * d* eigen';

rawX=mvnrnd(zeros(n,p),KK);
eta = rawX*beta;
Xmat = rawX * sqrt(6/var(eta)); % scale X to control the SNR (signal to noise ratio) to be 6 
Y = Xmat*beta + randn(n,1);

isNewSub = 50;       
isFPCA = 1;
[newXmat, stringed, disMatrix, x, t_x, SubId, xx, pred_x, pred_pc] = Stringing(Xmat, 'correlation', [], isNewSub, isFPCA, 1);
b = glmfit(pred_pc(SubId == 0), Y(SubId == 0), 'normal');
predy = glmval(b, pred_pc(SubId==1), 'identity');  
        
% compute error for the testing data
error = mean((predy-Y((ntrain+1):n)).^2);

