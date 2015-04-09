function [outtmp] = functmp(fnin, hin)
  
    indf = find(fnin~=0);
    indh = find(hin~=0);
    indout = intersect(indf, indh);
    outtmp1 = hin(indout).*log(hin(indout)./fnin(indout))-hin(indout)+fnin(indout);
    
%     ind2 = 1:length(outtmp1);
    indout1 = find(~isnan(outtmp1));
    indout2 = find(~isinf(outtmp1));
    outtmp = sum(outtmp1(intersect(indout1,indout2)));
    
      