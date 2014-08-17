function [invalid_data] =CheckData(y, t)
invalid_data = 0;
yy = cell2mat(y);
tt= cell2mat(t);
if (size(yy, 1) > 1) || (size(tt, 1) >1)
     fprintf(1, 'Error: the format of the data is incorrect!');
     invalid_data = 1;
     return;
 end
if (sum(isnan(yy))>0) || (sum(isnan(tt))>0)
      fprintf(1,'Error: FPCA is aborted because the data contain NaNs!');
      invalid_data = 1;
      return;
end
end

