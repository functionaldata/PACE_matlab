function [x_predOrig,y_predOrig,xi_est,yi_est]=getSOriCurves(x, y, t_x, t_y, mu_x, mu_y, sc_x, sc_y, lambda_x, lambda_y, lambda, sigma_x, sigma_y, sig1, nosvd, error, method, shrink, out1x, out1y, regular)

[muxSub, scxSub] = convertMuPhi(t_x, out1x, mu_x, sc_x, regular);
[muySub, scySub] = convertMuPhi(t_y, out1y, mu_y, sc_y, regular);
ncohort = length(y);
x_predOrig = cell(1, ncohort);
y_predOrig = cell(1, ncohort);

sigma1 = sig1;

if regular == 2 && strcmp(method,'CE')
    
    xx= reshape(cell2mat(x), length(x{1}), ncohort)';
    yy= reshape(cell2mat(y), length(y{1}), ncohort)';
    
    A22 = scxSub*lambda_x*scxSub'+error(1)*sigma1(1)*eye(length(x{1}));
    B22 = scxSub*diag(lambda)*scySub';
    D22 = scySub*lambda_y*scySub'+error(2)*sigma1(2)*eye(length(y{1}));
    A22i = pinv(A22-B22*pinv(D22)*B22');
    D22i = pinv(D22-B22'*pinv(A22)*B22);
    B22i = -pinv(A22)*B22*pinv(D22-B22'*pinv(A22)*B22);
    sig22i = [A22i B22i; B22i' D22i];
    sig12 = [lambda_x*scxSub' diag(lambda)*scySub'; diag(lambda)*scxSub' lambda_y*scySub'];
    %         sig11 = [lambda_x diag(lambda); diag(lambda) lambda_y];
    mu_ss = (sig12 * sig22i * ([xx yy] - repmat([muxSub muySub],ncohort,1))')';
    xi_est = mu_ss(:,1:nosvd);
    yi_est = mu_ss(:,(nosvd+1):end);
    
    x_predOrig = repmat(muxSub, ncohort,1)+xi_est*scxSub';
    x_predOrig = num2cell(x_predOrig,2)';
    y_predOrig = repmat(muySub, ncohort,1)+yi_est*scySub';
    y_predOrig = num2cell(y_predOrig,2)';
    
else
    xi_est = zeros(ncohort, nosvd);
    yi_est = zeros(ncohort, nosvd);
    zeta_xest = xi_est;
    zeta_yest = yi_est;
    scxi= scxSub;
    mux_i = muxSub;
    scyi= scySub;
    muy_i = muySub;
    for i = 1:ncohort
        
        if regular ~= 2
            scxi = scxSub{i};
            mux_i = muxSub{i};
            scyi = scySub{i};
            muy_i = muySub{i};
        end
        xi= x{i};
        yi= y{i};
        if strcmp(method,'CE')
            A22 = scxi*lambda_x*scxi'+error(1)*sigma1(1)*eye(length(xi));
            B22 = scxi*diag(lambda)*scyi';
            D22 = scyi*lambda_y*scyi'+error(2)*sigma1(2)*eye(length(yi));
            A22i = pinv(A22-B22*pinv(D22)*B22');
            D22i = pinv(D22-B22'*pinv(A22)*B22);
            B22i = -pinv(A22)*B22*pinv(D22-B22'*pinv(A22)*B22);
            sig22i = [A22i B22i; B22i' D22i];
            sig12 = [lambda_x*scxi' diag(lambda)*scyi'; diag(lambda)*scxi' lambda_y*scyi'];
            %             sig11 = [lambda_x diag(lambda); diag(lambda) lambda_y];
            mu_ss = sig12 * sig22i * ([xi yi] - [mux_i muy_i])';
            xi_est(i,:) = mu_ss(1:nosvd);
            yi_est(i,:) = mu_ss((nosvd+1):end);
            
        elseif strcmp(method,'IN')
            mx=length(xi);
            my=length(yi);
            for k=1:noeig
                prodx=(xi-mux_i).*scxi(:,k)';
                prody=(yi-muy_i).*scyi(:,k)';
                if error == 0
                    xi_est(i,k) = trapz(t_x{i},prodx);
                    yi_est(i,k) = trapz(t_y{i},prody);
                else
                    if shrink == 0
                        %xi_est(i,k) = romb(t{i},prod);
                        xi_est(i,k) = trapz(t_x{i},prodx);
                        yi_est(i,k) = trapz(t_y{i},prody);
                        %zeta_est(i,k)=trapzoid(prod,t{i});
                    else
                        %zeta_est(i,k) = romb(t{i},prod);
                        zeta_xest(i,k) = trapz(t_x{i},prodx);
                        zeta_yest(i,k) = trapz(t_y{i},prody);
                        xi_est(i,k)=lambda_x(k,k)*zeta_xest(i,k)/(lambda_x(k,k)+sigma_x/mx);
                        yi_est(i,k)=lambda_y(k,k)*zeta_yest(i,k)/(lambda_y(k,k)+sigma_y/my);
                    end
                end
            end
        end
        x_predOrig{i} = mux_i+xi_est(i,:)*scxi';
        y_predOrig{i} = muy_i+yi_est(i,:)*scyi';
    end
end
end
