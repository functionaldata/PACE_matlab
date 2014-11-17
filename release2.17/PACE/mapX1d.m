%function [newy] = mapX1d(x,y,newx)
%Map (x,y) to (newx,newy)
%x, y : vectors of 1 * n
%newx : vector of 1 * m
%newy : vector of 1*m
function [newy] = mapX1d(x,y,newx)

%      [ignore, loc] = ismember(newx, x);
     
     if (~all(ismember(newx, x)))
         warning('Interpolation occured: you might want to increase the out1 coverage')
     end
     if ((min(newx) + 100 * eps < min(x)) | (max(newx) > max(x) + 100 * eps))
         warning('Extrapolation occured')
     end
     newy = interp1(x, y, newx, 'linear', 'extrap');
     if (size(y, 2) == 1)
         newy = newy';
     end
     
%      if size(y,1) == 1
%         newy = y(loc);
%      elseif size(y,1) > 1
%         newy = y(loc,:);
%      else
%         fprintf(1,'y cannot be empty!\n');
%      end   
end
