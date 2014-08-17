function [res] = getRaw2dMean(tt, yy)
   tt = tt';
   out = unique(tt, 'rows');
   n = size(out, 1);
   count = zeros(1, n);
   y = zeros(1, n);
   for i = 1:n
       id = ismember(tt, out(i,:), 'rows');
       count(i) = sum(id);
       y(i) = sum(yy(id))/count(i);
   end
   tt = out';
   res = struct('tt',tt, 'yy', y, 'count',count);
end
