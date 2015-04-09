function [xhist_density,xhist_mids] = hist1D(x,nbins,h0,rangex)

    if ~isempty(nbins) && ~isempty(h0)
        fprintf(1,'Too many arguments! Only one of nbins and h0 is needed!\n');
        return;
    elseif ~isempty(nbins) && isempty(h0)
        h0 = range(rangex)/(nbins-1);
    elseif isempty(nbins) && ~isempty(h0)
        nbins = ceil(range(rangex)/h0+1);
        h0 = range(rangex)/(nbins-1);
    else
        fprintf(1,'Either nbins or h0 is needed here!\n');
        return;
    end
    
    histgrid = rangex(1)-h0/2:h0:rangex(2)+h0/2;
    xhist_n = histc(x,histgrid);
    xhist_density = xhist_n(1:end-1)/length(x)/h0;
    xhist_mids = rangex(1):h0:rangex(2);