function [ res ] = FPCAstep1(out1, noeig, ngrid, xcov)
out21 = linspace(min(out1),max(out1),ngrid);
h=range(out21)/(length(out21)-1);
opts.disp = 0;  %don't show the intermediate steps
ngrid = size(xcov,1);
[eigen d] = eigs(xcov,ngrid-2,'lm',opts);
d = d(~imag(d));
d = d(d > 0);
FVE = cumsum(d)/sum(d);
if isempty(noeig)
    noeig = length(d);
elseif noeig > length(d)
    noeig = length(d);
    fprintf(1,['Warning: at most ' num2str(noeig) ' number of PC can be selected!\n']);
end
eigen = eigen(:,1:noeig);
d = d(1:noeig);   

eigen = eigen./sqrt(h);
lambda = h*d';

for i = 1:noeig
    eigen(:,i) = eigen(:,i)/sqrt(trapz(out21,eigen(:,i).^2));    
    if eigen(2,i)< eigen(1,i)
       eigen(:,i)=-eigen(:,i);
    end
end
% interpolate from the normalized the eigenfunctions
phi=interp1(out21, eigen, out1, 'spline');
if noeig == 1
   phi = phi';
end
% normalize smoothed eigenfunctions
for i=1:noeig
   phi(:,i) = phi(:,i)/sqrt(trapz(out1,phi(:,i).^2));
end



Xnames = {'out1', 'phi', 'lambda', 'no_opt', 'FVE', 'names'};
res = {out1, phi, lambda,noeig, FVE, Xnames};


