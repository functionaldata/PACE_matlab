function[Y, group] = sim_obs(n, p1, lam1, lam0, phi1, phi0, PC_dist, sigma2)
% Simulate functional obserations with mean zero from a mixed distribution
%   n:          total number of curves generated
%   p1:         probability for generating population 1
%   lam1, lam0: row vectors of eigenvalues
%   phi1, phi0: matrices of eigenfunctions with the same size, with each
%               column standing for an eigenfunction
%   PC_dist:    the distribution of principal component scores. Support normal
%               and exponential
%   sigma2:     the variance of additional noise added to each observation 
%               point

n_grid = size(phi1, 1);
n1 = binornd(n, p1);
n0 = n - n1;
group = [ones(1, n1), zeros(1, n0)];
K = length(lam1);
if (strcmp(PC_dist, 'norm'))
    Xi = random(PC_dist, 0, 1, n, K);
elseif (strcmp(PC_dist, 'exp'))
    Xi = random(PC_dist, 1, n, K) - 1; % mean 0 and var 1
end

Xi(1:n1, :) = Xi(1:n1, :) * diag(sqrt(lam1));
Xi((n1 + 1):n, :) = Xi((n1 + 1):n, :) * diag(sqrt(lam0));

Y = nan * ones(n, n_grid);
Y(1:n1, :) = Xi(1:n1, :) * phi1';
Y((n1 + 1):n, :) = Xi((n1 + 1):n, :) * phi0';

% Add noise
Y = Y + sqrt(sigma2) .* randn(n, n_grid);

% randomize the order
ord = randsample(n, n);
Y = Y(ord, :);
group = group(ord);

end