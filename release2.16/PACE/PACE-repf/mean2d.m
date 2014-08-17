function [invalid,mu]=mean2d(bw,kernel,xin,yin,win, out, count)
% use local weighted least squares kernel smoothing.
% input bw(1,2):    specified bandwidth, b(1) for x axis b(2) for y axis
% input kernel:     specified kernel function, candidates:'epan','rectangle'
% input xin(2,n):   matrix of predictors. xin(:,i) is the coordinate of the ith
%                   predictor.
% input yin(n,1):   vector of responses. y(i) is the response valued at x(:,i)
% input win(1,n):   vector of case weights. win(i) is the weight of the ith 
%                   observation.
%  input count: the number of repeats for each unique point in xin.
% deleting observations with weights 0
active=find(win~=0);
xin=xin(:,active);
yin=yin(active);
win=win(active);
invalid=0;
%gap=[];
mu = zeros(1, size(out, 2));
for i=1:size(out, 2)
   xnow =out(:,i);
   [invalid temp] = GetEstimate(xin, yin, win, bw, kernel, xnow, count); 
   if invalid == 1
        mu = [];
        return;
   end
   mu(i) = temp;
end
  