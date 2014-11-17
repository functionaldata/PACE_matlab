% manifold kernel smoothing, predict functional values/FPC scores based on
% functional manifold components
%
% y - m-by-d matrix of functional manifold components
% S - a structure, see maniMDS
% kernel - smoothing kernel [default epan]
% h - smoothing bandwidth [default by Cross-Validation]
% out - leave-one-out indicator
% d - manifold dimensions used [default size(y,2)]
% score - indicator, use 1 if predict FPC scores
% x - predicted functional values/FPC scores
% hind - optional indicator - 0: perform bandwidth and Knn search when h empty [default]
%                             1: perform only bandwidth search when h empty
%                             Note: for faster computing, specify a value for h;
%                               in this case, no warning results when y is not
%                               sufficiently close to the data manifold scores S.y.  

function x = maniKS(y,S,kernel,h,out,d,score,hind)

if nargin<8|isempty(hind) hind = 0; end
if nargin<7|isempty(score) score=0; end
if nargin<6|isempty(d) d=size(y,2); end
if nargin<5|isempty(out) out=0; end
if nargin<4 h=[]; end
if nargin<3|isempty(kernel) kernel='epan'; end

% distances between points to be reverse-mapped (y) and sample points (S.Y)
d_mani = L2_distance(y(:,1:d)',S.Y(:,1:d)',0);

if isempty(h)
	[knnind,knnD] = knnsearch(S.Y,y,'K',3);
	if any(any(S.h < knnD))
		warning('Warn:maniKS',...
		  'Some input manifold components may be too far from the data;\n curves obtained from maniKS() may be distorted')
	end
	hmin = max(max(knnD)); % Ensures at least three data points are within one bandwidth for all components
	h = bwCV(S,kernel,hmin,hind);      
end


if isa(h,'char') % For K-nearest neighbor search
    K = min(str2num(h),S.N-out);
    if out == 1
        outnum = sum(d_mani'==0);
    else
        outnum = zeros(1,size(y,1));
    end
    [Sort_dis,Ord] = sort(d_mani');
    if K==1 
        x = S.X_reg(Ord(1,:)+outnum,:);
        return; % simply return closest trajectory
    end
    w = zeros(size(y, 1), S.N);
	for l = 1:size(w, 1)
		idx = (outnum(l) + 1):(K + outnum(l)); % locations of K-nearest neighbors, excluding identical points
		w(l, Ord(idx, l)) = 1./Sort_dis(idx, l); % weight the neighbors by inverse distance
	end		
else % bandwidth search instead of Knn
    H = repmat(h,[1,size(y,1)]);
	w = kernelval(repmat(1./H',[1,S.N]).*d_mani,kernel);
end

if out==1 w(find(w==kernelval(0,kernel)))=0; end
if score==1
    X = S.xi;
else
    X = S.X_reg;
end
if all(sum(w,2)) % check for zero weights
    x = repmat(1./sum(w,2),[1,S.N]).*w*X;
else
    idx = find(sum(w,2)>0); % nonzero weights
    x = repmat(NaN,[size(y,1),S.M]);
    x(idx,:) = repmat(1./sum(w(idx,:),2),[1,S.N]).*w(idx,:)*X; % smooth sample trajectories
end

end