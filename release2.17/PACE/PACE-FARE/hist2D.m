function [xhist2D,hout,histtmp] = hist2D(x2D, nbins, gridsize,xrange)
    
    if ~isempty(nbins) && ~isempty(gridsize)
        fprintf(1,'Too many arguments! Only one of nbins and gridsize is needed!\n');
        return;
    elseif ~isempty(nbins) && isempty(gridsize)
        delta = range(xrange)/(nbins-1);
    elseif isempty(nbins) && ~isempty(gridsize)
        nbins = ceil(range(xrange)/gridsize+1);
        delta = range(xrange)/(nbins-1);
    else
        fprintf(1,'Either nbins or histgrid is needed here!\n');
        return;
    end
    
    hout = xrange(1):delta:xrange(2);
    count = hist3(x2D',{hout,hout});
    count = count+count';
    histtmp = count/(size(x2D,2)*delta^2);
%     xhist2D = histtmp/(sum(histtmp(:))*hx*hy);
    xhist2D = histtmp/trapz2(histtmp,hout,hout);
     