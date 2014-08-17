% ===========
% Description:
% ===========
%             This is the main function to perform dynamical correlation.
%             Reference: Dubin,J.,Muller,H. (2005), "Dynamical Correlation 
%             for Multivariate Longitudinal Data".
%
%======
%Usage:
%======
%
% function [dyncorr]=dynCorr(y,t,regular)
%
%================
%Input Arguments:
%================
%
%Input y:         m*n cell array, y{i,j} is the measurements for the ith 
%                 function(curve) and jth subject, i=1,...,m, j=1,...,n.
%Input t:         m*n cell array, t{i,j} is the vector of time points for
%                 the ith function(curve) and jth subject for which corresponding
%                 measurements y{i,j} are available.
%Input regular:   0, sparse (or irregular) functional data.      
%                 1, regular data with missing values.
%                 2, completely balanced (regular) data.
%
%================
%Output Arguments:
%================
%
%Output dyncorr:  m*m matrix, dyncorr(i,j) is the dynamical correlation
%                 between ith and jth function(curve).



function [dyncorr]=dynCorr(y,t,regular)

if regular~=0 && regular~=1 && regular ~=2
     fprintf(1,'Error: Please specify design (regular = )');
     dyncorr=NaN;
     return;


elseif (regular==0) || (regular==1)
    
    p=setOptions('regular',regular,'kernel','epan');
    ndep=size(y,1);
    nsub=size(y,2);
    newy=cell(ndep,nsub);
    for i=1:ndep
        
        [yy] = FPCA(y(i,:),t(i,:),p);
        newy(i,:) = yy{17};
    end
    y=newy;
end

[dyncorr] = dynCorr1(y);

for i=2:ndep
    for j=1:(i-1)
        dyncorr(i,j)=dyncorr(j,i);
    end
end


        
    
