function [res] = getRawRepfCov(y, t, muforcov, out1, out2)
[tx1 tx2 tx3] = meshgrid(out2, out2, out1);
tx1 = tx1(:)';
tx2 = tx2(:)';
tx3 = tx3(:)';
cyy = zeros(1, length(tx1));
count = cyy;
n = length(y);
for i = 1:n
   for j = 1:length(out1)
       [tf, ids] = ismember(t{i}(1,:), out1(j));
       if sum(ids)>0 
          ynow = y{i}(ids==1)-muforcov(:,j)';     
          [y1, y2] = meshgrid(ynow);
          y1 = y1(:)';
          y2 = y2(:)';
          idin = (tx3 == out1(j));
          cyyidin = y1.*y2;
          cyy(idin) = cyy(idin) + cyyidin;
          count(idin) = count(idin)+1;
       end
   end
end
  tpairn = [tx1; tx2; tx3];
  tpairn = tpairn(:, count > 0);
  cyy =cyy(count>0);
  count = count(count>0);
  cyy = cyy./count;
  res = struct('tpairn',tpairn, 'cyy',cyy,'count',count);
end
  