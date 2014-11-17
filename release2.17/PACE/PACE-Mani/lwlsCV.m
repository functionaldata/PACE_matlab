% local weighted linear smoothing function, used for interpolation for
% missing values and noise removing
%
% X,T - cell array with observed data and time points
% X_reg,t - estimated trajectories on balanced time grid t
% kernel - smoothing kernel [default epan]
% h - kernel bandwidth [default by cross-validation]

function [X_reg,t] = lwlsCV(X,T,t,kernel,h)

if nargin<5 h = []; end
if nargin<4|isempty(kernel) kernel = 'epan'; end
if nargin<3 t = unique(cell2mat(T)); end

N = length(X);
for k=1:N
    tmptX=X{k};
    tmptT=T{k};
    if isempty(h)
        r = range(tmptT);
        n = length(tmptT);
        h0 = max(tmptT(2:n)-tmptT(1:(n-1)))*1.01;
        if strcmp(kernel,'guass') h0 = h0/2; end
        q = (r/(4*h0))^(1/9);
        bw = sort(q.^(0:9).*h0);
    else
        bw = h;
    end
    gcv = Inf*ones(1,length(bw));
    for i=1:length(bw)
        mu = lwlval(tmptX,tmptT,tmptT(2:(n-1)),kernel,bw(i),1);
        gcv(i) = (tmptX(2:(n-1))-mu)*(tmptX(2:(n-1))-mu)';
        clear mu;
    end
    [val,ind] = min(gcv);
    X_reg(k,:) = lwlval(tmptX,tmptT,t,kernel,bw(ind),0);
    clear tmptT tmptX bw;
end

end

