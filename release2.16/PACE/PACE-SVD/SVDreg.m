%function [b newy] = SVDreg(res, newx, new_tx, y, new_ty, isYFun, K)
%======
%Input:
%======
%  xx, yy: The returned values from FPCreg. See FPCreg() for more details.
%  newx:   1*numNewSub cell array for predictor function x,  newx{i} is the row vector of 
%          measurements for the ith new subject, i=1,...,numNewSub.
%  new_tx: 1*numNewSub cell array, new_tx{i} is the row vector of time points for the ith 
%          new subject on which corresponding measurements newx{i} are taken, 
%          i=1,...,numNewSub.
%  y:      the scalar response Y vector used for estimation 
%          Only needed when Y is a scalar.
%  new_ty: 1*numNewSub cell array, new_ty{i} is the row vector of time points
%          for the ith new subject on which corresponding measurements newy{i} will be taken, 
%          i = 1,...,numNewSub, or [] if y is scalar response.
%  K:      positive integer; number of principal components of predictor x
%          used in prediction must be smaller than or equal to that used in
%          regression.
%=======
%Output:
%=======
%  newy:   1*numNewSub cell array for response y,  newy{i} is the response
%          for the i-th new subject, i = 1,..., numNewSub.

function [b newy] = SVDreg(res, newx, new_tx, y, new_ty, isYFun, K)


    xx=getVal(res,'xx');
    yy=getVal(res,'yy');
    ss_x = getVal(res, 'ss_x');
    ss_y = getVal(res, 'ss_y');
%     rho = getVal(res, 'rho_opt');
    
    if isempty(K)
        K=size(ss_x,2);
    end
    
    b = getB(ss_x, ss_y, K, K, isYFun);
    out1 = getVal(xx, 'out1copy');
    out_x = getVal(res, 'out_x');
    sigma = getVal(xx, 'rho_opt');
    sc_x = getVal(res, 'sc_x');
    lambda_x = getVal(res,'lambda_x')';
    lambda_x = lambda_x(1:K,1:K);
    nsub = length(newx);
            
    if getVal(xx,'regular') == 2
        new_tx = new_tx{1};
        newxmu = interp1(out1, getVal(xx, 'mucopy'), new_tx, 'spline');
        if K==1
            newxsc = interp1(out_x, sc_x, new_tx, 'spline')';
        else
            newxsc = interp1(out_x, sc_x, new_tx, 'spline');
        end
        newxcov = newxsc * lambda_x * newxsc' + sigma * eye(length(new_tx));
        newxscore = (lambda_x * newxsc'* pinv(newxcov) * (reshape(cell2mat(newx),length(new_tx), nsub)-repmat(newxmu,nsub,1)'))';
    else
        newxscore = ones(nsub, K);
        for j=1:nsub
            newtxi = new_tx{j};
            newxmu = interp1(out1, getVal(xx, 'mucopy'), newtxi, 'spline');
            if K==1
                newxsc = interp1(out_x, sc_x, newtxi, 'spline')';
            else
                newxsc = interp1(out_x, sc_x, newtxi, 'spline');
            end
            newxcov = newxsc * (lambda_x) * newxsc' + sigma * eye(length(newtxi));
            newxscore(j,:) = (lambda_x * newxsc'* pinv(newxcov) * (newx{j}-newxmu)')';
        end
    end
    
    if isYFun == 1
        b = b(1:K,1:K);
        sc_y = getVal(res, 'sc_y');
        out_y = getVal(res, 'out_y');
        if getVal(yy,'regular') == 2
            new_ty = new_ty{1};
            newymu = interp1(getVal(yy, 'out1copy'), getVal(yy, 'mucopy'), new_ty, 'spline');
            if K==1
                newysc = interp1(out_y, sc_y, new_ty, 'spline')';
            else
                newysc = interp1(out_y, sc_y, new_ty, 'spline');
            end
            newy = mat2cell(repmat(newymu,nsub,1) + newxscore * b * newysc',ones(1,nsub))';
        else
            newy = cell(1,nsub);
            for j = 1:nsub
                newtyi = new_ty{j};
                newymu = interp1(getVal(yy, 'out1copy'), getVal(yy, 'mucopy'), newtyi, 'spline');
                sc_y = getVal(res, 'sc_y');
                if K==1
                    newysc = interp1(out_y, sc_y, newtyi, 'spline')';
                else
                    newysc = interp1(out_y, sc_y, newtyi, 'spline');
                end
                newy{j} = newymu + newxscore(j,:) * b * newysc';
            end
        end
    else
        b = b(1:K);
        newy = mat2cell(mean(y)*ones(nsub,1) + newxscore * b',ones(1,nsub))';
    end
end
