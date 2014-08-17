
function [X_pred,t,mse1,mse2] = predX(S,kernel,h,out,d,X0)

if nargin<6 X0 = []; end
if nargin<5|isempty(d) d = size(S.Y,2); end
if nargin<4|isempty(out) out = 0; end
if nargin<2|isempty(kernel) kernel = 'epan'; end
if nargin<3|isempty(h) h = h10cv(S,kernel); end
    
t = S.t;
X_pred = repmat(NaN,[length(S.Outliers),S.M]);
X_pred(~S.Outliers,:) = maniKS(S.Y,S,kernel,h,out,d);
if isempty(X0)
    mse1 = NaN;
    mse2 = NaN;
elseif nargout>2
    idx = ~isnan(X_pred(:,1));
    idx2 = ~isnan(X_pred(~S.Outliers,1));
    if sum(idx)
        mse1 = trace(L2_distance(X_pred(idx,:)',X0(idx,:)',0).^2*range(S.t)/(S.M-1))/sum(idx);
        mse2 = trace(L2_distance(S.X_reg(idx2,:)',X0(idx,:)',0).^2*range(S.t)/(S.M-1))/sum(idx2);
    else
        mse1 = NaN;
        mse2 = NaN;
    end
end
        
end