% function [fx,outgrid,bw,xhist_density] = hist2fx(x, nbins, h0, bw, kern, outgrid)
% hist2fx is a function that first generates the histogram of the data and then
% smoothes the histogram to get a density estimate fx

function [fx,bw,xhist_density,xhist_mids] = hist2fx(x, nbins, h0, bw, kern, outgrid,rangex)

    % set default values for truncate & outgrid
    if isempty(kern)
        kern = 'gauss';
    end
    
    if isempty(rangex)
        rangex = [min(x),max(x)];
    end
    
    [xhist_density,xhist_mids] = hist1D(x,nbins,h0,rangex);
    
    if isempty(outgrid)
        outgrid = xhist_mids;
    end
        
    if bw == 0
        bw = gcv_lwls(xhist_density,xhist_mids,kern,1,1,0,0,'off');
    end
    
    [invalid, fxtmp] = lwls(bw,kern,1,1,0,xhist_mids,xhist_density',ones(1,length(xhist_mids)),outgrid);
    
    ind = find(~isinf(fxtmp) & ~isnan(fxtmp));
    fxtmp = interp1(outgrid(ind),fxtmp(ind),outgrid,'spline');
    fxtmp(fxtmp < 0) = 0;
    fx = fxtmp/trapz(outgrid,fxtmp);
    
     
