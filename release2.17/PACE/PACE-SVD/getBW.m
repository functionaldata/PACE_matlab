function [res invalid]=getBW(x,t_x,xx, y,t_y,yy,rcov,tpair, regular)

if nargin<=4
    regular=0;
end
% error=1;

%initialize
nsvd = 3;

% n = length(x);
% p = setOptions('regular',regular,'selection_k', nsvd, 'screePlot',0);
% xx = FPCA(x,t_x,p);
% p = setOptions('regular',regular,'selection_k', nsvd, 'screePlot',0);
% yy = FPCA(y,t_y,p);
ngrid=51;
% sigma1_x = getVal(xx,'rho_opt');
% sigma1_y = getVal(yy,'rho_opt');
%obtain smoothed cross covariance
% x_mu=getVal(xx,'mucopy');
% y_mu=getVal(yy,'mucopy');
% out_x=getVal(xx,'out1copy');
% out_y=getVal(yy,'out1copy');
out2x=getVal(xx,'out21copy');
out2y=getVal(yy,'out21copy');
% [rcov, tpair] = getRawCCOV(x,t_x,y,t_y,x_mu,y_mu,1,regular);
win=ones(1,length(tpair));

bwgcv = gcv2_mullwlsn(t_x,t_y,ngrid,regular,'gauss',tpair,rcov,win);
bwcddt = [bwgcv(1)*(0.6:0.2:3)' bwgcv(2)*(0.6:0.2:3)'];
nbw = length(bwcddt);
ccov_s = cell(1,nbw);
sc_x = cell(1,nbw);
sc_y = cell(1,nbw);
ccovfit = cell(1,nbw);
lambda = cell(2,nbw);
lambda_x = cell(2,nbw);
lambda_y = cell(2,nbw);
invalid = zeros(1, nbw);
for ibw =1:nbw
    [invalid(ibw), ccov]=mullwlsk_2(bwcddt(ibw,:),'gauss',tpair,rcov',win,out2x,out2y);
    if invalid(ibw) == 0
        ccov_s{ibw} = ccov;
        [sc_xi di sc_yi] = svds(ccov, nsvd);
%         save temp;
        hx = trapz(out2x, sc_xi.^2)';
        hy = trapz(out2y,  sc_yi.^2)';
        sc_xi = sc_xi / diag(sqrt(hx));
        sc_yi = sc_yi / diag(sqrt(hy));
        di = diag(di).*sqrt(hx.*hy);
        di = di(1:nsvd);
        lambda{1,ibw} = di;
        lambda_xi = (sc_xi' * getVal(xx,'xcovfitcopy') * sc_xi);
        lambda_yi = (sc_yi' * getVal(yy,'xcovfitcopy') * sc_yi);
        lambda_xi = lambda_xi .* (hx*hx');
        lambda_yi = lambda_yi .* (hy*hy');
        lambda_x{1,ibw} = diag(lambda_xi(1:nsvd,1:nsvd));
        lambda_y{1,ibw} = diag(lambda_yi(1:nsvd,1:nsvd));
%         lambda_xi = zeros(1,nsvd);
%         lambda_yi = zeros(1,nsvd);
%         for isvd = 1:nsvd
%             z= repmat(sc_xi(:,isvd),1,length(out2x)) .* getVal(xx,'xcovfitcopy') .* repmat(sc_xi(:,isvd)',length(out2x),1);
%             lambda_xi(isvd)=trapz2(z, out2x, out2x);
%             z= repmat(sc_yi(:,isvd),1,length(out2y)) .* getVal(yy,'xcovfitcopy') .* repmat(sc_yi(:,isvd)',length(out2y),1);
%             lambda_yi(isvd)=trapz2(z, out2y, out2y);
%         end
%         lambda_x{2,ibw} = lambda_xi(1:nsvd);
%         lambda_y{2,ibw} = lambda_yi(1:nsvd);
%         sc_xi = sc_xi(:,1:nsvd);
%         sc_yi = sc_yi(:,1:nsvd);
%         sc_x{ibw} = sc_xi;
%         sc_y{ibw} = sc_yi;
%         ccovfit{ibw} = sc_xi * diag(di) * sc_yi';
%         lambdai = zeros(1,nsvd);
%         for isvd = 1:nsvd
%            z= repmat(sc_xi(:,isvd),1,size(ccov_s{ibw},2)) .* ccov_s{ibw} .* repmat(sc_yi(:,isvd)',size(ccov_s{ibw},1),1); 
%            lambdai(isvd)=trapz2(z, out2x, out2y);
%         end
%         lambda{1,ibw} = lambdai';
    end
end

tmp = 1:nbw;
validbw = tmp(invalid ==0);

for ibw =validbw
    sc_xi = sc_x{ibw};
    sc_yi = sc_y{ibw};
    for jbw = validbw
        lambdai = zeros(1,nsvd);
        for isvd = 1:nsvd
           z= repmat(sc_xi(:,isvd),1,size(ccov_s{jbw},2)) .* ccov_s{jbw} .* repmat(sc_yi(:,isvd)',size(ccov_s{jbw},1),1); 
           lambdai(isvd)=trapz2(z, out2x, out2y);
        end
        lambda{jbw,ibw} = lambdai';
    end
end

lambda2 = cell(2,nbw);
for ibw =validbw
    sc_xi = sc_x{ibw};
    sc_yi = sc_y{ibw};
    sc_xi = interp1(out2x, sc_xi, out_x, 'spline');
    sc_yi = interp1(out2y, sc_yi, out_y, 'spline');
    for i=1:size(sc_xi,2)
        sc_xi(:,i) = sc_xi(:,i)/sqrt(trapz(out_x,sc_xi(:,i).^2));
    end
    for i=1:size(sc_yi,2)
        sc_yi(:,i) = sc_yi(:,i)/sqrt(trapz(out_y,sc_yi(:,i).^2));
    end
    for j = 1:2
        lambdai = lambda{j,ibw};
        ccovfiti = sc_xi * diag(lambdai) * sc_yi';
        res1names = {'sigma1_x' 'sigma1_y' 'xx' 'yy' 'sc_x', 'sc_y', 'out_x', 'out_y', 'out2x', 'out2y', 'ccovfit', 'lambda', 'lambda_x', 'lambda_y', 'names'};
        res1 = {sigma1_x, sigma1_y, xx, yy, sc_xi, sc_yi, getVal(xx, 'out1copy'), getVal(yy, 'out1copy'), out2x, out2y, ccovfiti, lambdai, lambda_x{2,ibw}, lambda_y{2,ibw}, res1names};
        [xi_est, yi_est, ss_var]=getSScores(res1, x, t_x, y, t_y, nsvd, [1 1], 'CE', 1,  regular, 'cv');
        cov_cond = zeros(n,nsvd);
        for i=1:n
            for k=1:nsvd
                cov_cond(i,k)=ss_var{i}(k,nsvd+k);
            end
        end
        ss_cov = cov([xi_est yi_est]);
        sv=zeros(1,nsvd);
        for k=1:nsvd
            sv(k)=mean(cov_cond(:,k))+ss_cov(k,nsvd+k);
        end
        lambda2{j,ibw} = sv;
    end
end

resnames = {'sc_x', 'sc_y', 'out2x', 'out2y' 'ccov_s' 'ccovfit', 'lambda',  'lambda_x', 'lambda_y', 'bwgcv', 'bwcddt', 'rcov', 'tpair','names'};
res = {sc_x, sc_y, out2x, out2y, ccov_s, ccovfit, lambda, lambda_x, lambda_y,  bwgcv, bwcddt, rcov, tpair, resnames};

% cov_cond = zeros(n,nsvd);
% for i=1:n
%     for k=1:nsvd
%         cov_cond(i,k)=ss_var{i}(k,nsvd+k);
%     end
% end
% ss_cov = cov([xi_est yi_est]);
% sv=[];
% for k=1:nsvd
%     sv(k)=mean(cov_cond(:,k))+ss_cov(k,nsvd+k);
% end
% fc = sv(1)/sqrt(lambda_x(1)*lambda_y(1));
% sv_mean = mean(cov_cond(:,1));
% sv_cov = ss_cov(1,nsvd+1);


% [ss_x]=getScores(x, t_x, getVal(xx,'mu'), sc_x, lambda_x, getVal(xx,'sigma'), nsvd, error, 'CE', 0, getVal(xx,'out1copy'), regular, 'cv');
% [ss_y]=getScores(y, t_y, getVal(yy,'mu'), sc_y, lambda_y, getVal(yy,'sigma'), nsvd, error, 'CE', 0, getVal(yy,'out1copy'), regular, 'cv');

% resnames = ['fc' 'sv' 'sv_mean' 'sv_cov' 'ss_x' 'ss_y' 'xi_est' 'yi_est' 'ss_var' 'rho_opt' 'sig1' getVal(res,'names')];
% res = [fc sv sv_mean sv_cov ss_x ss_y xi_est yi_est {ss_var} rho_opt sig1 res];
% res{end} = resnames;
%ctime=cputime()-ctime;
%fprintf(1,['Final Choice: \n' num2str(ec) '=[' num2str(bwccov(1)) ',' num2str(bwccov(2)) '],' num2str(selex) ',' num2str(seley) '@' num2str(ctime) 's\n']);
end
