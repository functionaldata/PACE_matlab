%% [mu,dyn_grid]=nonlindyn(t,y,method,tout,N_f)
% This function computes the empirical dynamics of (t,y(t),y'(t))
% Input t: 1 * N cell of observed time points. Need to be coherent for all subjects
% Input y: 1 * N cell of observed repeated measurements for each subject, corresponding to each cell of t.
% Input method: a character string for the kernel to be used, default is 'gauss'
%               'epan' : epanechikov kernel
%               'gauss' : gaussian kernel
%               'gausvar' : variant of gaussian kernel
%               'rect' : rectangular kernel
%               'quar' : quartic kernel
% Input tout: m * 1 vector of the output time grid ; default is equadistance 1*100 grid between min(t{1}) and max(t{1}).               
% Input N_f:  scalar that gives the output size of trajectory y grid; default is 200.
% Output mu:  m * N_f matirx of the estimated dynamic surface on given time grid and trajectory grid (N_f equidistant grid on the support region)
% Output dyn_grid: m * N_f matrix of the output grid corresponding to mu.
   
function [mu,dyn_grid] = nonlindyn(t,y,method,tout,N_f)
    if nargin==4
     N_f=200;
   end
    if nargin == 3 
    stp=range(t{1})/100;
    tout=min(t{1}):stp:max(t{1});
    N_f=200;
    end            
  
     if nargin==2
     stp=range(t{1})/100;
    tout=min(t{1}):stp:max(t{1});
    N_f=200;          
    method='gauss';
    end
   
    tout_len=length(tout);
    n=length(y);
    
    % change non-cell object t into cells.
    if(~iscell(t))
        tmp=t;
        t=cell(1,n);
        for i=1:n
            t{i}=tmp;
        end
    end
    
    
    y_CellSize=size(y{1});
    if(y_CellSize(1)==1)
       y=y';
    end
    
% obtaining the optimal bandwidth for x0 and x1 smoothing.      
     r = range(t{1});     
     N = length(t{1});  %number of all observed time points
     h0 = minb(t{1},2)*1.5; %minimum candidate bandwidth for x0 smoothing
     h1=minb(t{1},3)*1.5;     % minimum candidate bandwidth for x1 smoothing
     h0 = min(h0, r);
     h1=min(h1,r);
     
      q0 = (r/(4*h0))^(1/9);
       q1 = (r/(4*h1))^(1/9);
     [out1,ignore, idx] = unique(t{1});   %distinct  observed time points
     %fprintf(1,'Bandwidth choices for mean function: \n');
     bw0 = sort(q0.^(0:9).*h0) ;                %create 10 h candidates
     bw1 = sort(q1.^(0:9).*h1) ;                %create 10 h candidates
     k0 = mykernel(0, method);
     t_mat=kron(t{1},ones(1,n));    
     y_mat=cell2mat(y);
     yy=y_mat(:);


t_i=t{1};   % input t for each subject
t1_i=(t{1}(2:end)+t{1}(1:end-1))/2; %input t for raw derivative
gcv_sum0=zeros(1,length(bw0));
gcv_sum1=zeros(1,length(bw1));
for k=1:length(bw0)
  for i=1:n
     y_i=(y_mat(i,:))';
     w=ones(1,length(t_i))/bw1(k);
     w0=ones(1,length(t_i))/bw0(k);
    [invalid,mu]=lwls(bw0(k),method,1,1,0,t_i,y_i,w0,t_i,0);
    y1=y_i(2:end)-y_i(1:end-1);
    [invalid,mu1]=lwls(bw1(k),method,2,2,1,t_i,y_i,w,t1_i',0);
    gcv_sum0(k)=gcv_sum0(k)+(y_i'-mu)*(y_i'-mu)';
    gcv_sum1(k)=gcv_sum1(k)+(y1'-mu1)*(y1'-mu1)';
  end
  gcv_sum0(k)=gcv_sum0(k)/bw0(k);
  gcv_sum1(k)=gcv_sum1(k)/bw1(k);
end

 clear h0 h1;
 h0=bw0(gcv_sum0==min(gcv_sum0));
 h1=bw1(gcv_sum1==min(gcv_sum1));
 
 
 %  Smooth estimates for x0 and x1
x_0=zeros(n,tout_len);% will be used to store  smoothed trajectories
x_1=zeros(n,tout_len);% will be used to store smoothed derivatives
 for i=1:n
     y_i=(y_mat(i,:));
     w=ones(1,length(t_i))/h1;
     w0=ones(1,length(t_i))/h0;
     [invalid,mu]=lwls(h0,method,1,1,0,t_i,y_i',w0,tout,0);
     [invalid,mu1]=lwls(h1,method,2,2,1,t_i,y_i',w,tout,0);
    x_0(i,:)=mu;
     x_1(i,:)=mu1;
  end
  
   
% find the support for estimating the surface
mmin_x0=min(x_0(:));
mmax_x0=max(x_0(:));
mmin_x0=ceil(mmin_x0*100)/100;
mmax_x0=floor(mmax_x0*100)/100;
f=zeros(N_f,2,tout_len);  % N_f is the number of grid for x(t), tout_len is the length of output time grid.

Z=zeros(n,2,tout_len);  % Z  is the residual matrix
 for j=1:tout_len
   f(:,1,j)=linspace(mmin_x0,mmax_x0,N_f);  %  first dimension is the grid of x_pred
 end
gridf=linspace(mmin_x0,mmax_x0,N_f);  % grid of X(t) where f might be estimated

a=size(x_0);
a=a(1)*a(2);
b=size(x_1);
b=b(1)*b(2);
XX=reshape(x_0',1,a);         
YY=reshape(x_1',1,b);
[bopt2,gcv]=gcv_lwls(YY,XX,method,1,1,0,2,'on',0); 
  
 
bw=zeros(tout_len,1);  %bandwidth used at each time point

% estimating the dynamic surface f
f_x=NaN(n,tout_len);
 I_x0=NaN(n,tout_len);
 for k=1:tout_len
    xk=(x_0(:,k))';
    x1k=x_1(:,k);
    x11k=x1k;
    posMinxk=find(xk==min(xk));
    xk(posMinxk)=[];
    x11k(posMinxk)=[];
    posMaxxk=find(xk==max(xk));
    xk(posMaxxk)=[];
    x11k(posMaxxk)=[];
    ngridf=find(gridf<=max(xk) & gridf>=min(xk));  
    xgridf=gridf(ngridf);
    if size(x11k,1)~=1
       x11k=x11k';
    end
    bopt2temp=gcv_lwls(x11k,xk,method,1,1,0,0,'off',0); % at different time points the bandwidth for estimating f varies.
    bw(k)=bopt2temp*1.1;  % or choose bandwidth at each time point?
    w=ones(1,n)/bopt2temp;
    [invalid,mu2]=lwls(bopt2temp,method,1,1,0,(x_0(:,k))',x1k,w,xgridf,0); 
    [x0k_sorted I_x0k]=sort(x_0(:,k));
    I_x0(:,k)=I_x0k;
    [invalid2,f2]=lwls(bopt2temp,method,1,1,0,(x_0(:,k))',x1k,w,x0k_sorted,0);  
     f_x(:,k)=f2;    %% f_x: each column stores the est. of f w.r.t observed x (after smoothing).
                %% in reality, given each time t, the est. of f w.r.t biggest/smallest x_0(.,t) does not 
                %% have a reasonable estimate, thus to be removed in estimation of sos_h,sos_x.
     f(ngridf,2,k)=mu2;
end
    mu=f(:,2,:);
    mu=squeeze(mu);
    dyn_grid=squeeze(f(:,1,:));
    


