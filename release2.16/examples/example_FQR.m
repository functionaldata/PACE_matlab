%% Time
% 100 subjects
base=(0.02:0.02:0.98)'; % length 49
T=repmat(base,[100,1]); % length 49*100
Ni=49*ones(1,100);
Tcell=mat2cell(T',1,Ni); % 1*100 cell

%% Covariate
zbase=randn(100,2);
z=repmat(reshape(zbase,[],1),[1,49]); % 2-dim independent covariates, normal(0,1)
z=reshape(z',[],2); % 4900*2

%% Principal decomposition

err=randn(4900,2)*0.05; % measurement errors for principal scores
% principal scores with eigenvalue 1 and 0.09.
A1=(sqrt(2)/2)*z(:,1)+(sqrt(2)/2)*z(:,2)+err(:,1); % first index: sqrt(2)/2*z1+sqrt(2)/2*z2; alpha_1(x)=x;
A2=0.3*((sqrt(2)/2)*z(:,1)-(sqrt(2)/2)*z(:,2))+err(:,2); % second index: sqrt(2)/2*z1-sqrt(2)/2*z2; alpha_2(x)=0.3*x;

% eigenfunctions
phi1=-sqrt(2)*cos(pi*(T+0.5)); % vector of length 4900
phi2=sqrt(2)*sin(pi*(T+0.5)); % same


%% Response (no measurement errors)
X=1+(T-0.5).^2+A1.*phi1+A2.*phi2;
  
Xcell=mat2cell(X',1,Ni);

%% Fit

p = path;
isExist = regexp(p, 'PACE');
if isempty(isExist) == 1
  addpath(genpath('../PACE/'));
end
p1 = setOptions('yname','x', 'selection_k', 2, 'verbose','off');  


[S,W]=FQR( Xcell, Tcell, zbase, p1, 1);

%% Compare with true values
%Note here that eigenvectors and eigenvalues are unique up to a sign
%change.


r=-2:0.01:3;
figure(1);
hold on;
plot(r,r,'color','green');
legend('Observations','Estimated link function','True link function');

t=-3:0.01:4;
figure(3);
hold on;
plot(t,-0.3*t,'color','green');
legend('Observations','Estimated link function','True link function');


