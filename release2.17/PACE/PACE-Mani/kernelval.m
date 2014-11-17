% function calculating the kernel values
%
% x - time locations, vector or matrix
% kernel: a chacracter string
% 'quar' - quartic kernel
% 'epan' - epanechikov kernel
% 'gauss' - gaussian kernel
% 'gausvar' - variant of gaussian kernel
% 'rect' - rectangular kernel
% kx: kernel values

function kx = kernelval(x,kernel)

if strcmp(kernel,'quar')
    kx = (15/16).*(abs(x) <= 1).*(1-x.*x).^2;
elseif strcmp(kernel,'epan')
    kx = 0.75.*(abs(x) <= 1).*(1-x.*x);
elseif strcmp(kernel,'rect')
    kx = 0.5.*(abs(x) <= 1);
elseif strcmp(kernel,'gausvar')
    kx = (1/sqrt(2*pi)).*exp(-0.5.*x.*x).*(1.25-0.25*x.^2);
else
    kx = (1/sqrt(2*pi)).*exp(-0.5.*x.*x);
end

end