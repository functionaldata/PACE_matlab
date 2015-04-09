% function [f2D,out1,out2,denstmp] = hist2fx2D(x2D, nbins, gridsize, bw, kern, out1, out2)
% hist2fx2D is a function that first generates the 2-D histogram of the data and then
% smoothes the histogram to get a 2-D density estimate f2D

function [f2D,bw2,denstmp] = hist2fx2D(x,x2D,nbins,gridsize,bw2,kern,ngrid1,out21,xrange)

    if isempty(kern)
        kern = 'gauss';
    end
    
    [denstmp,mids,histtmp] = hist2D(x2D,nbins,gridsize,xrange);

    yin = histtmp(:);
    [xin1,xin2] = meshgrid(mids,mids);
    xin = [xin2(:),xin1(:)]';
    win = ones(1,length(yin));
    x2dtmp = struct('tpairn',xin,'cxxn',yin','indx',[],'win',win,'cyy',yin','count',[]);
        
    if sum(bw2) == 0
        bw2 = gcv_mullwlsn(num2cell(mids),ngrid1,2,1,kern,x2dtmp,'off');
    end
    [invalid, f2Dtmp] = mullwlsk(bw2,kern,xin,yin,win,out21,out21);
    
    f2Dtmp = (f2Dtmp+f2Dtmp')/2;
    f2Dtmp(f2Dtmp<0) = 0;
    f2D = f2Dtmp/trapz2(f2Dtmp,out21,out21);
    
    