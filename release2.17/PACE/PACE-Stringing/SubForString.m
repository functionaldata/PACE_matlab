function [SubId, invalid] = SubForString(Xmat, isNewSub)
   n = size(Xmat, 1);
   invalid = 0; 
   SubId = ones(1,n);
   if isempty(isNewSub) || all(isNewSub == 0)
      SubId = zeros(1,n);
   elseif length(isNewSub) == 1 && isNewSub > 0            %when isNewSub is scalar, the last number 
      if isNewSub == n
         fprintf(1,'Error: invalid isNewSub, at least one or more subjects are not for predicton!\n');
         invalid=1;
         return;
      end 
      SubId(1:(n-isNewSub)) =0 ;
   elseif length(isNewSub) == n          %when isNewSub is an indicator,1 : new subject , 0: observed subject
     if all(isNewSub == 1)
        fprintf(1,'Error: invalid isNewSub, at least one or more subjects are not for predicton!\n');
        invalid=1;
        return;
     end
     SubId(isNewSub == 0) = 0;
   else
      fprintf(1,'Warning: isNewSub must be a positive integer or a 0-1 indicator vector!Reset isNewSub = [] now!\n');
   end