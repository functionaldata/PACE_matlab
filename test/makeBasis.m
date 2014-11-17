function [ basis_true1, basis_true0 ] = makeBasis( grid_pts, K_true, shift)
%makeBasis: Make a (shifted) Fourier basis of length K_true supported on grid_pts
% Assume the basis functions are supported on [0, 1]
% Input: grid_pts need to be between [0, 1]
% Output: basis_true0 is shifted while basis_true1 is not

n_grid = length(grid_pts);
basis_true1 = zeros(n_grid, K_true);
basis_true1(:, 1) = 1;
basis_true0 = basis_true1;
for (i = 2:K_true) 
    if (mod(i, 2) == 0)
        basis_true1(:, i) = sqrt(2) * cos(2 * pi * i / 2 * grid_pts);
        basis_true0(:, i) = sqrt(2) * cos(2 * pi * i / 2 * (grid_pts - shift));
    else
        basis_true1(:, i) = sqrt(2) * sin(2 * pi * (i - 1) / 2 * grid_pts);
        basis_true0(:, i) = sqrt(2) * sin(2 * pi * (i - 1) / 2 * (grid_pts - shift));
    end
end  

end

