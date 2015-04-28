%  Example for DiffEq.m

%  Random design, 200 subjects with 1 to 8 measurements from each subject.
%  Time interval is [0,1]. The number of measurements for each subject is
%  uniformly distributed on {1,2,...,8}. The timepoints for each subject
%  are distributed as Beta(0.4,0.3).


%generate data set
%clear all;

p = path;
isExist = regexp(p, 'PACE');
if isempty(isExist) == 1
  addpath(genpath('../PACE/'));
end
rng(123, 'twister');
mtp = 8;        %at most 8 repeated measurements in the simulated data
ncohort=200;    %200 subjects in the simulated data
lint=1;
y=cell(1,ncohort);
t=cell(1,ncohort);
newy = y;
newt = t;
xi=zeros(ncohort,2);
ngrid = 100;

%Case iii) regular data with missing values (regular = 1)   
%ni_missing = zeros(1,ncohort);
%for i = 1:ncohort
%  ni_missing(i) = poissrnd(mysample(1:mtp,1,0),1,1); 
%  if ni_missing(i) >= mtp
%     ni_missing(i) = mtp-1;
% end
%end

for i=1:ncohort

   %Case i) Sparse and irregular case (regular = 0) 
   ntp=ceil(mtp*rand(1));
   t{i}=lint*rand(1,ntp);                 
   %newt{i} = lint*rand(1,ntp);
   newt{i} = sort(betarnd(0.4,0.3,1,ntp));

   %Case ii) complete balance case (regular = 2) 
   %t{i} = linspace(0,lint,mtp);                      
   %newt{i} = lint*rand(1,ntp);

   %Case iii) regular data with missing values (regular = 1)   
   %t{i} = linspace(0,lint,mtp);                    
   %newt{i} = linspace(0,lint,mtp);
   %if ni_missing(i) > 0
   %  id_missing = mysample(1:mtp, ni_missing(i),0);
   %  t{i}(id_missing) = [];
   %  newt{i}(mysample(1:mtp,ni_missing(i),0)) = [];
   %end

   xi(i,:)=[3*randn(1) 2*randn(1)];     %generate 2 Principal components
                                        %1st PC score: N(0,9), 2nd PC score: N(0,4)
   %generate the repeated measurements with measurement errors
   
   y{i}=genMeanFun_1(t{i},0)+xi(i,:)*genEigenFun_1(t{i},0)+randn(1,length(t{i}));
   newy{i} = genMeanFun_1(newt{i},0) + xi(i,:)*genEigenFun_1(newt{i},0)+randn(1,length(newt{i}));
          
   %measurement error is distributed as N(0,1)
end


 
%========================================================================================================
%Run FPCAder
p = setDerOptions('yname','x','selection_k', 'FVE','FVE_threshold', 0.85,'screePlot',0, 'corrPlot',0,...
		  'designPlot',0,'numBins',0, 'nder',0:1,'newdata',linspace(0,1,ngrid), 'ngrid',ngrid,'verbose','on');  

% Here, the options "p" is set as the follows: yname = 'x', use FVE method for 
% selection of No. of FPCs (default is FVE, but if only the curve itself is of interest, set it as 'BIC1'), 
% FVE_threshold = 0.85, 
% create scree plot and design plot (by default, they are not displayed),
% do not perform prebinning for the input data (by default, the program try to bin the data)
% input data, perform curve estimation as well as 1st and 2nd derivative estimation 
% (by default nder = 0, only the curve itself will be estimated),
% define output time grids as 100 equidistant grids in [0,1],
% use 100 grids for the covariance surface and their partial derivative of 
% covariance surface, display diagnostic messages,
% and the rest uses default values                           
% It is important that 'nder' starts with 0, meaning the estimation of the trajectories themselves.
% By default, 'verbose' is set to be 'on'. You can suppress the diagnostic messages by setting it to be 'off'.

%Use FPCder() to recover the functional object for y and the results is a cell array
%that can be assessible by getVal(), e.g., to get eigenfunctions "phi", use
%phi = getVal(yy, 'phi');

fprintf(1,'Recover the individual trajectories using FPCder():\n');
time=cputime;
[XXd] = FPCder(y,t,p);
time = cputime-time;
display(['Total time in FPCder fitting is : ' num2str(time) ' sec.']);

%==========================================================================
%Run DiffEq.m

[toutx beta1 zdc Gz R2 phiz VarXder VarZ] = DiffEq(XXd);
n=ncohort;

%plot all the Z(i)
figure(1)
for i=1:n
    plot(toutx,zdc(i,:),'k-');
    hold on;
end
xlabel('Time','fontsize',14);
ylabel('Z_i(t)','fontsize',14)

%plot \beta(t) and eigenfunctions of Z(t)
figure(2)
subplot(1,2,1)
plot(toutx,beta1,'k-')
del=range(beta1)/20;

xlabel('Time','fontsize',14);
ylabel('Deynamic transfer function \beta(t)','fontsize',14);
subplot(1,2,2)
plot(toutx,phiz(:,1),'k-',toutx,phiz(:,2),'k--');
xlabel('Time','fontsize',14);
ylabel('Eigenfunctions of Z(t)','fontsize',14);

%plot Variance Functions
figure(3)
subplot(1,2,1)
plot(toutx,VarXder,'k--',toutx,VarZ,'k-');
xlabel('Time','fontsize',14)
ylabel('Variance functions of X^{(1)}(t) and Z(t)','fontsize',14)
subplot(1,2,2)
plot(toutx,R2,'k-')
xlabel('Time','fontsize',14)
ylabel('R^2(t)','fontsize',14)
