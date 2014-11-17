% function beta = transferXY(xx, yy, nder, gridin, gridout)
% Calculates transfer functions and cross-relationships between Gaussian
% processes X and Y. \beta(s,t) = E(Y^{(\nu)}(t)|X(s))
% Reference: M\"uller, H.G. and Yang, W. (2009). Dynamic relations for
% sparsely sampled Gaussian processes. TEST.
%======
%Input:
%======
%  xx, yy:  the returned values from FPCder. See FPCder() for more details.
%  nder:    the order \nu derivative of Y at t that is regressed on X(s).
%  gridin:  optional, the input grid of the transfer function, corresponds to argument s.
%  gridout: optional, the output grid of the transfer function, corresponds to argument t.
%
%=======
%Output:
%=======
%  beta:    transfer function, a matrix whose dimension is determined by
%           gridin and gridout.
function [beta Sigma] = transferXY(xx, yy, nder, gridin, gridout)
    if nargin < 5
        if nargin == 3
            gridin = getVal(xx,'out21copy');
        else
            fprintf(1,'Error: Not enough input arguments!');
            beta = [];
            return;
        end
        gridout = gridin;
    end
    x = getVal(xx,'y');
    t_x = getVal(xx,'t');
    y = getVal(yy,'y');
    t_y = getVal(yy,'t');
    ngrid = 51;
    regular = min(getVal(xx,'regular'),getVal(yy,'regular'));
    lambda = getVal(xx,'lambda');
    phi_x = getVal(xx,'eigencopy');
    if iscell(phi_x)
        phi_x = phi_x{1};
    end
    phi_y = getVal(yy,'eigencopy');
    if ~iscell(phi_y) && nder == 0
        phi0_y = phi_y;
    elseif ~iscell(phi_y) && nder ~= 0
        fprintf(1,'Error: yy must be the output from FPCder()!');
        beta = [];
        return;
    else
        phi0_y = phi_y{1};
        phi_y = phi_y{nder+1};
    end
    x_mu=getVal(xx,'mucopy');
    if iscell(x_mu)
        x_mu=x_mu{1};
    end    
    y_mu=getVal(yy,'mucopy');
    if iscell(y_mu)
        y_mu=y_mu{1};
    end
    out2x=getVal(xx,'out21copy');
    out2y=getVal(yy,'out21copy');
    [rccov, tpair] = getRawCCOV(x,t_x,y,t_y,x_mu,y_mu,1,regular);
    win=ones(1,length(tpair));
    bwccov = gcv2_mullwlsn(t_x,t_y,ngrid,regular,'gauss',tpair,rccov,win);
    [invalid, ccov_s]=mullwlsk_2(bwccov,'gauss',tpair,rccov',win,out2x,out2y);

    Sigma = (phi_x'*ccov_s*phi0_y);
    phis = interp1(out2x,phi_x, gridin, 'spline');
    phit = interp1(out2y,phi_y, gridout,'spline');
    beta = (phis*Sigma*phit') ./repmat((phis.^2)*lambda',1,length(gridout));
end
