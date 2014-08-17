% function to find local weighted linear fits
%
% x,t - observed data and time vectors
% y0,t0 - estimated trajectory y0 on time grid t0
% kernel - smoothing kernel
% h - kernel bandwidth
% out - leave-one-out indicator

function y0 = lwlval(x,t,t0,kernel,h,out)

D = repmat(t0',[1,length(t)])-repmat(t,[length(t0),1]);
W = kernelval(D/h,kernel);
if out==1 W(find(W==kernelval(0,kernel)))=0; end
for k=1:length(t0)
    X2 = [ones(1,length(t));t-t0(k)];
    tmpt = pinv(X2*diag(W(k,:))*X2')*X2*diag(W(k,:))*x';
    y0(k) = tmpt(1);
    clear X2 tmpt;
end

end