function [lambda, phi, covfit, eigen] = geteigen(xcov,out1,out21,varargin)

h = range(out21)/(length(out21)-1);
opts.disp = 0;  %don't show the intermediate steps
ngrid = size(xcov,1);

[eigen d] = eigs(xcov,ngrid-2,'lm',opts);

d = diag(d);
iddx = find(imag(d));                           %remove any imaginary eigenvalues                                               
if ~isempty(iddx) 
   error('Error: %d eigenvalues are complex. The estimated auto-covariance surface is not symmetric!',length(iddx));
end                                                                                                                      
iddx = find(d <= 0);                                                                                                       
if ~isempty(iddx)   
   if length(varargin) > 0                                                                                                     
     fprintf(1,['Warning: ' num2str(length(iddx)) ' real eigenvalues are negative or zero and are removed!\n']);
   end                    
   d = d(d > 0);                                  %retain only the positive eigenvalues
   eigen = eigen(:,d > 0);                                  
end                                              

noeig = size(eigen,2);
eigen = eigen./sqrt(h);
lambda = h*d';

% normalize smoothed eigenfunctions
for i=1:noeig
   eigen(:,i) = eigen(:,i)/sqrt(trapz(out21,eigen(:,i).^2));
   if eigen(2,i) < eigen(1,i)
      eigen(:,i) = -eigen(:,i);
    end
end

phi=interp1(out21', eigen, out1', 'spline');
for i=1:noeig
   phi(:,i) = phi(:,i)/sqrt(trapz(out1,phi(:,i).^2));
end

covfit = eigen*diag(lambda)*eigen';



