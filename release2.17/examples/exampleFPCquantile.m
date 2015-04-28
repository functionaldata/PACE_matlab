% Example for FPCquantile.m 
%
% Random design, 200 subjects, 30 measurements on each subject.
% 150 for estimation and 50 for prediction.
% time domain is [0,10] for predictor X and the observations for X are taken equidistantly.
% Y conditioning on X has a Gaussian distribution. 

p = path;
isExist = regexp(p, 'PACE');
if isempty(isExist) == 1
addpath(genpath('../PACE/'));
end

n = 200;
npred = 50;
rng(123, 'twister');

numK = 4; % true number of components used to generate data
x =cell(1,n);
t_x =cell(1,n);    
xi =zeros(n, numK);
lint=10;
y = 1:n;
c1 = 0.5;
sqrtlambda = [4, 3, 2.75, 2.25];
for i = 1:n
    mtp = 30;
    t_x{i} = linspace(0,lint,mtp);
    for j = 1:numK
        xi(i,j)= sqrtlambda(j)*randn(1);  
    end 
    x{i}= mu_true(t_x{i},lint)+xi(i,:)*xeig(t_x{i},lint, numK)+ c1*randn(1,length(t_x{i}));
    y(i)= 5*randn(1)+ 2*sum(xi(i,:));
end

outQ =  [0.1, 0.25, 0.5, 0.75, 0.9, 0.95];
param_X = setOptions('regular',2,'verbose','on');
isNewSub = npred;
time=cputime;
[predQ, predF, xx] = FPCquantile(x, t_x, y, outQ, param_X, isNewSub);
time = cputime-time;
display(['Total time in FPCquantile is : ' num2str(time) ' sec.']);

% Calculate the true conditional quantiles for the testing data.
T = zeros(npred, length(outQ));
testing = (n-npred+1) : n;
for i= 1:npred
    T(i,:)= icdf('Normal',outQ, 2*sum(xi(testing(i), :)), 5);       
end

%calculate the mean squared error for the testing data
abserror = zeros(length(outQ), 1);
for j = 1:length(outQ)
    abserror(j) = mean(abs(T(:,j)- predQ(testing, j)));
end
%plot the true against the predicted conditional quantiles for 
%all testing subjects.
k=0;
for i = 1:6
    k= k+1;
    subplot(3,2, k)
    plot(T(:,i), predQ(testing,i), '.',  [min(T(:,i)), max(T(:,i))], [min(T(:,i)), max(T(:,i))])
    xlabel(['q(', num2str(outQ(i)), ')'] , 'Interpreter','latex', 'FontSize', 12)
    ylabel(['$$\hat q($$', num2str(outQ(i)), ')'] , 'Interpreter','latex', 'FontSize', 12)
end
  
