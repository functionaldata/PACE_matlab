% function [res] = fpcsdaic(x,bw1,bw2,covmethod,kern,xrange,selection_k,...
%                  FVE_threshold,nbins,nbins1,ngrid,ngrid1,newdata,control,screePlot)
% Function fpcsd stands for Funtional Principal Component (FPC) Analysis for
% Sparse Density Estimation f(x)
%
% Input:
%     x           1*n cell, sparse observations for each subject.
%     bw1         bandwidth to be used in estimating the overall mean density function,
%                 bw1 = 0 will use the GCV bandwidth selection.
%     bw2         bandwidth for the 2-dimensional density function f2d(x,y) (when covmethod=0)
%                 or the covariance matrix Ghat(x,y)(when covmethod = 1). 
%                 bw2 = c(0, 0) will use the GCV bandwidth selection.
%   covmethod     method to be used in estimating the covariance surface Ghat(x,y)
%                 0  covariance calculated directly from the smoothed 2-D
%                    density function f2d(x,y) and mean function fmuhat(x):
%                    Ghat(x,y) = f2d(x,y)-fmuhat(x)*fmuhat(y).
%                 1  construct the histgram estimates f2d_hist and fmu_hist
%                    for 2-D density function and mean function
%                    respectively and then obtain smoothed covariance
%                    Ghat(x,y) by applying 2-D smoothing to
%                    f2d_hist(x,y)-fmu_hist(x)*fmu_hist(y). 
%     kern        a character string to define the kernel to be used in the
%                 1-D or 2-D smoothing
%                 kernel = 'epan'  ==> Epanechnikov kernel 
%                          'rect'  ==> Rectangular kernel
%                          'gauss'  ==> Gaussian kernel   [Default]
%     xrange      the range of the output x, a vector containing the minimum and the maximum.
%                 if [] then the range of the input x will be used.
%  selection_k    the method of choosing the number of principal components K.
%                 'PPIC': use pseudo-Poisson information criterion.   [Default]
%                 'FVE' (fraction of variance explained) : use scree plot 
%                         approach to select number of principal components),                  
%                 see "FVE_threshold" below
%                 positive integer K: user-specified number of principal
%                 components.
% FVE_threshold   a positive number that is between 0 and 1 [Default is 0.85]
%                 It is used with the option selection_k = 'FVE' to select
%                 the number of principal components that explain at least
%                 "FVE_threshold" of total variation (the fraction
%                 of variance explained).
%     nbins       number of bins to be used to construct the 1-D histogram of pooled data.
%                 [Default is 101]
%     nbins1      number of bins to be used to construct the individual histograms.
%                 [Default is 101]
%     ngrid       number of support points in each direction of covariance surface.      
%                 [Default is 51]               
%     ngrid1      number of support points in the GCV procudere (selecting bandwidth) for 
%                 the 2-D density (when covmethod = 0) or covariance
%                 surface (when covmethod = 1).      [Default is 30]
%     newdata     a row vector of user-defined output time grids for all curves. 
%                 If newdata = [], then "out1" corresponds to the set of distinct
%                 time points from the pooled data.
%                 "newdata" is supposed to be a vector in ascending order on
%                  the domain of the functions.         [Default is []]
%     control     'auto', Select K by minimizing PPIC score, or find the 
%                         first K that exceeds the FVE_threshold.  [Default]
%                 'look', a scree plot (FVE% Vs No. of PC) will be generated based 
%                         on K <= 15. User will be prompted to enter user-specified
%                         K after viewing scree plot. The feature can be combined 
%                         with any setting of selection_k.
%   screePlot     indicator of whether to create the scree plot
%                 1  a scree plot will be created         
%                 0  no scree plot will be created      [Default]
% 
% Output: A cell of estimates
%     fhat        N*nn matrix, estimated subject-specific density functions valued at distinct
%                 input time points from all subjects with ascending order, corresponding to out1.
%     fmuhat      1*N vector, estimated mean density function valued at distinct input time 
%                 points from all subjects with ascending order, corresponding to out1.
%     Ghat        ngrid*ngrid matrix, smoothed covariance surface using the method "covmethod", 
%                 corresponding to out21.
%     covfit      ngrid*ngrid matrix, fitted covariance surface, based on truncated estimate
%                 of eigenvalues and eigenfunctions, corresponding to out21.              
%     f2d         ngrid*ngrid matrix, smoothed 2-D density function when covmethod = 0, 
%                 corresponding to out21.
%     ifhat       N*nn matrix, estimated subject-specific intensity functions, corresponding to out1.
%     tauhat      1*nn vector, number of observations for each subject.
%     out1        1*N vector, distinct input time points from all subjects with ascending order 
%                 if newdata = []; otherwise, it is the same as newdata.
%     out21       1*ngrid vector, a grid of time points for which the smoothed covariance surface
%                 assumes values.
%     noeig       integer, automatically or subjectively selected value of K, the number of selected components. 
%     lambda      estimates of the eigenvalues.
%     phi         estimates of the eigenfunctions, corresponding to out1.
%     xihat       nn*K matrix, estimates of the FPC scores.
%     bw          a structure that contains bandwidths bw1, bw2 used in the estimation and 
%                 the indicator covmethod.
%     id          1*nn vector, the index of subjects with at least two observations.
%     loglike     log-likelihood obtained when choosing different values for K, the number of components 
%                 (available when selection_k = "PPIC").
%     fhist       N*nn matrix, subject-specific density functions estimated by smoothing individual 
%                 histograms, corresponding to out1 (available when selection_k = 'PPIC'). 

function [res] = fpcsdaic(x,bw1,bw2,covmethod,kern,xrange,selection_k,FVE_threshold,nbins,nbins1,ngrid,ngrid1,newdata,control,screePlot)

    xx = cell2mat(x);
    % set default values
    if isempty(bw1)
        bw1 = 0;
    end
    if isempty(bw2)
        bw2 = [0,0];
    end
    if isempty(covmethod)
        covmethod = 1;
    end
    if isempty(kern)
        kern = 'gauss';
    end
    if isempty(xrange)
        xrange = [min(xx),max(xx)];
    end
    if isempty(selection_k)
        selection_k = 'PPIC';
    end
    if strcmp(selection_k,'FVE') == 1 && isempty(FVE_threshold)
        FVE_threshold = 0.85;
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
    if isempty(control)
        control = 'auto';
    end
    if isempty(screePlot)
        screePlot = 0;
    end
        
    % Truncate the data if necessary
    nn1 = length(x);
    if min(xx) < xrange(1) || max(xx) > xrange(2)
        for i = 1:nn1
            x{i} = x{i}(x{i}<=xrange(2) & x{i}>=xrange(1));
        end
    end
        
    % Require at least 2 samples per subject, if not satisfied, then eliminate the subjects
    lengthtmp = zeros(1,nn1);
    for i = 1:nn1
        lengthtmp(i) = length(x{i});
    end
    ind_geq2 = find(lengthtmp>=2);
    tauhat = lengthtmp(ind_geq2);
    x = x(ind_geq2);
    nn = length(x);
    if nn < nn1
        fprintf(1,'NOT all data are used!\n');
    end
    
    % Estimate the mean density function f_mu and generate the bandwidth bw1
    xx = cell2mat(x);
    if isempty(newdata)
        out1 = sort(unique(xx));
    else
        out1 = newdata;
    end
    [f_mu,bw1] = hist2fx(xx,nbins,[],bw1,kern,out1,xrange); 
      
    % Estimate the covariance surface G and eigenvalues and eigenvectors
    out21 = xrange(1):diff(xrange)/(ngrid-1):xrange(2);
    
    x2Dtmp = cell(1,nn);
    for i = 1:nn
        x2Dtmp{i} = nchoosek(x{i},2)';
    end
    x2D = cell2mat(x2Dtmp);
    
    if covmethod == 0
        [f2D,bw2] = hist2fx2D(x,x2D,ngrid,[],bw2,kern,ngrid1,out21,xrange);
        fmu21 = interp1(out1,f_mu,out21,'spline');
        Ghat = f2D-fmu21'*fmu21;
    else
        [mutmp,outtmp] = hist1D(xx,ngrid,[],xrange);
        [f2dtmp] = hist2D(x2D,ngrid,[],xrange);
        Ghattmp = f2dtmp-mutmp'*mutmp;
        yin3 = Ghattmp(:);
        [xin31,xin32] = meshgrid(outtmp);
        xin3 = [xin32(:),xin31(:)]';
        clear xin31 xin32;
        indneq = find(xin3(1,:)~=xin3(2,:));
        yin3 = yin3(indneq);
        xin3 = xin3(:,indneq);
        win3 = ones(1,length(yin3));
        Ghattmp = struct('tpairn',xin3,'cxxn',yin3','indx',[],'win',win3,'cyy',yin3','count',[]);
        if sum(bw2) == 0
            bw2 = gcv_mullwlsn(num2cell(outtmp),ngrid1,2,1,kern,Ghattmp,'off');        
        end
        [invalid,Ghat] = mullwlsk(bw2,kern,xin3,yin3,win3,out21,out21);
        f2D = [];
    end
    Ghat = (Ghat+Ghat')/2;
    
    [lambda,phi,covfit] = geteigen(Ghat,out1,out21,1);
    noeigtmp = size(phi,2);    

    % Estimate the FPC scores xi and the density function f for each subject
    % The columns in fhat correspond to the density function estimates
    % Calculate the unadjusted density   
    
    FVE = cumsum(lambda)/sum(lambda);
    noeigCopy= find(FVE > FVE_threshold,1,'first');
    pc_options = {'PPIC','FVE','user'};

    if ischar(selection_k)
        
        k_id = strmatch(selection_k, pc_options,'exact');
        if isempty(k_id)
            fprintf(1,['Warning: Invalid method name for selection_k! Reset to "FVE" method with threshold = ' num2str(FVE_threshold) '\n']);
            k_id = 2;
        end
        if k_id == 1
            [fhat,xihat,noeig,loglike,fhist] = getfhat(x,xrange,f_mu,phi,out1,0,noeigtmp,kern,nbins1);
        elseif k_id == 2
            noeig = noeigCopy;
            [fhat,xihat,noeig,loglike,fhist] = getfhat(x,xrange,f_mu,phi,out1,noeig,noeigtmp,kern,nbins1);
        end
        
    elseif isnumeric(selection_k) && selection_k > 0
        
        if selection_k > noeigtmp;
            noeig = noeigtmp;
            fprintf(1,['Warning: at most ' num2str(noeigtmp) ' number of PC can be selected!\n']);
        else
            noeig = selection_k;
        end
        k_id = 3;
        [fhat,xihat,noeig,loglike,fhist] = getfhat(x,xrange,f_mu,phi,out1,noeig,noeigtmp,kern,nbins1);
        
    else
        
        fprintf(1,['Error: "selection_k" must be a positive integer! Reset to "FVE" method with threshold = ' num2str(FVE_threshold) '\n']);
        k_id = 2;
        noeig = noeigCopy;
        [fhat,xihat,noeig,loglike,fhist] = getfhat(x,xrange,f_mu,phi,out1,noeig,noeigtmp,kern,nbins1);
        
    end    
    
    fprintf(1, ['Best number of principal components selected by ' pc_options{k_id} ': ' num2str(noeig) '.\n']) ;
    if k_id ~= 2
         fprintf(1,['It accounts for ' num2str(roundoff(FVE(noeig), 4)*100) '%% of total variation.\n']);
    else
         fprintf(1,['It accounts for ' num2str(roundoff(FVE(noeig), 4)*100) '%% of total variation (threshold = ' num2str(FVE_threshold) ').\n']);
    end
    fprintf(1,['FVE calculated from ' num2str(noeigtmp) ' possible eigenvalues: \n']);
    disp(FVE);
    
    
    if strcmp(control, 'look')
       fve2plot = FVE(1:ceil(noeigtmp/2))*100;
       figure;
       plot(0:length(fve2plot), [0 fve2plot],'--ro', 'LineWidth',2,...
       'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'g', 'MarkerSize',5);
       xlabel('\bf{No. of Principal Components}');
       ylabel('\bf{FVE (%)}');
       title(['\bf{Fraction of variance explained by No. of PC (threshold = ' num2str(FVE_threshold) ')}'])
       hold on
       plot(linspace(0,noeigCopy,30),ones(30,1)*FVE(noeigCopy)*100,'b',...
       ones(30,1)*noeigCopy, linspace(0,FVE(noeigCopy)*100,30),'b');
       text(noeigCopy+0.2, FVE(noeigCopy)*100-10, {[' k = ' num2str(noeigCopy) ', FVE = ' num2str(roundoff(FVE(noeigCopy)*100,3)) '%'] ' (threshold choice)'});
       axis([0 length(FVE)+1 0 101]);
       hold off

       noeig = input('Enter the number of principal components you want to choose:\nK=');
       fprintf(1, ['You just chose ' num2str(noeig) ' principal component(s).\n']);
       fprintf(1,['It accounts for ' num2str(roundoff(FVE(noeig), 4)*100) '%% of total variation.\n\n']);
    end

    
    %Now output scree plot based on final no_opt
    if screePlot == 1
      figure;
      plot(0:length(FVE), [0 FVE]*100, '--ro', 'LineWidth',2,'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'g', 'MarkerSize',5);
      xlabel('\bf{No. of Principal Components}');
      ylabel('\bf{FVE (%)}');
      title('\bf{Fraction of variance explained by No. of PC}')
      hold on;
      plot(linspace(0,noeig,30), ones(30,1)*FVE(noeig)*100,'b',...
      ones(30,1)*noeig, linspace(0,FVE(noeig)*100,30),'b');
      text(noeig+0.2, FVE(noeig)*100-10, {[' k = ' num2str(noeig) ', FVE = ' num2str(roundoff(FVE(noeig)*100,3)) '%'] ' (final choice)'});
      axis([0 length(FVE)+1 0 101]);
      hold off;
    end
    
    % Calculate the intensity function
    ifhat = fhat*diag(tauhat);
    
            
    % Define the outputs
    bw = struct('bw1',bw1,'bw2',bw2,'covmethod',covmethod);
        
    resNames = {'fhat','fmuhat','Ghat','covfit','f2d','ifhat','tauhat','out1',...
        'out21','noeig','lambda','phi','xihat','bw','id','loglike','fhist'};
    res = {fhat,f_mu,Ghat,covfit,f2D,ifhat,tauhat,out1,out21,noeig,lambda,phi,...
        xihat,bw,ind_geq2,loglike,fhist,resNames};




