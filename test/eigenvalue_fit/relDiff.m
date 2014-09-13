function [reldiff] = relDiff(xin, xtrue)
% Returns the relative difference of xin and xtrue.
% xtrue is a vector and xin can be a matrix. Each row of xin is compared
% against the true value xtrue.

if (size(xtrue, 1) > size(xtrue, 2))
    xtrue = xtrue';
end

n = size(xin, 1);
reldiff = (xin - ones(n, 1) * xtrue) ./ (ones(n, 1) * xtrue);

end