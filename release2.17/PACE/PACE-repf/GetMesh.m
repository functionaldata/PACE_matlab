function [ out ] = GetMesh(x1, x2)
[y1 y2] = meshgrid(x1, x2);
y1 = reshape(y1, 1, size(y1,1)*size(y1,2));
y2 = reshape(y2, 1, size(y2,1)*size(y2,2));
out = [y1;y2];
end

