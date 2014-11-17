% 10 fold cross validation function to choose bandwidth h
%
% S - a structure, see maniMDS
% kernel - smoothing kernel
% hmin - optional minimum bandwidth
% hind - optional indicator - 0: perform bandwidth and Knn search, 1: perform only bandwidth search
% bw - selected bandwidth

function bw = bwCV(S,kernel,hmin,hind)

if nargin<4|isempty(hind) hind = 0; end
if nargin<3 hmin = []; end

[h,Hcv] = h10cv(S,kernel,hmin);
if hind
	bw = h;
else
	[K,Kcv] = hknn10cv(S,kernel);
	if min(Kcv{2})<min(Hcv{2})
		bw = K;
	else
		bw = h;
	end
end


end          
        