function [xi_est, yi_est, ss_var, rho_opt, sig1, x_predOrig, y_predOrig, cv]=getSScores(res, x, t_x, y, t_y, nsvd, error, method, shrink,  regular, rho)

    mu_x=getVal(getVal(res,'xx'),'mucopy');
    sc_x=getVal(res,'sc_x');
    lambda_x = getVal(res,'lambda_x');
    sigma_x = getVal(getVal(res,'xx'),'sigma');
    out1x = getVal(res,'out_x');
    mu_y=getVal(getVal(res,'yy'),'mucopy');
    sc_y=getVal(res,'sc_y');
    lambda_y = getVal(res,'lambda_y');
    sigma_y = getVal(getVal(res,'yy'),'sigma');
    out1y = getVal(res,'out_y');
    rho_opt=[];        

    n = length(x);
    sig1 = [sigma_x sigma_y];

    %Compute iterative residuals
    if rho ~= -1
      for j = 1:2
        SE_i = inf*ones(n,2);
        %get fitted curves based on t{i} for a given TOL(j)
        [xo yo]=getSOriCurves(res, x, y, t_x, t_y, mu_x, mu_y, sc_x, sc_y, lambda_x, lambda_y, getVal(res,'lambda'),sigma_x, sigma_y, sig1, nsvd, error, method, shrink, out1x, out1y, regular);
%         xo = getOriCurves(x, t_x, mu_x, sc_x, lambda_x, sigma_x, sig1(1), nsvd, error(1), method, shrink, out1x, regular);        
%         yo = getOriCurves(y, t_y, mu_y, sc_y, lambda_y, sigma_y, sig1(2), nsvd, error(2), method, shrink, out1y, regular);        
        for i = 1:n
           SE_i(i,1) = mean((x{i}-xo{i}).^2);
           SE_i(i,2) = mean((y{i}-yo{i}).^2);
        end
        sig1 = mean(SE_i);
      end
    end

    if isnumeric(rho) 
       if rho >= 0 || rho == -1
         rho_opt = rho;
       else
	 fprintf(1,'Warning: rho should not be negative! Reset it to cv choice now!\n');
         rho = 'cv';
       end
    end

    if ischar(rho) && strcmp('cv',rho)

        %Compute gamma
        %fprintf(1, ['The new sigma after 2 iterations is ' num2str(sig1) '\n']);
        T_x = range(cell2mat(t_x));
        gamma_x = ((trapz(out1x, mu_x.^2)+sum(diag(lambda_x)))/T_x)^(0.5);
        T_y = range(cell2mat(t_y));
        gamma_y = ((trapz(out1y, mu_y.^2)+sum(diag(lambda_y)))/T_y)^(0.5);
        %fprintf(1,['Gamma = ' num2str(gamma) '\n']);
        %alpha = linspace(0.001*gamma, 0.05*gamma, 50);
        %alpha = linspace(0.005,0.22,50);
        rhomax = [gamma_x gamma_y];
        rhomin = sig1;
        rho = [];
        rho(:,1) = linspace(rhomin(1),rhomax(1),50);
        rho(:,2) = linspace(rhomin(2),rhomax(2),50);

        ni = zeros(n,2);
        tjID = zeros(n,2);
        for i = 1:n
            for k=1:2
                ni(i,1) = length(t_x{i});
                ni(i,2) = length(t_y{i});
                if ni(i,k) > 1
                    tjID(i,k) = mysample(1:ni(i,k),1,0); %find random index of random time point to be used in the CV prediction part
                    %for ni >= 2
                end
            end
        end

        %find optimal rho from subjects with ni >= 2
        if all(ni(:,1) == 1) || any(error == 0) || all(ni(:,2) == 1) 
           rho_opt = min(rho);
           cv = [];
        else
           [rho_opt cv] = cv_srho(res, x, t_x, y, t_y, nsvd, sig1, method, shrink, regular, rho, ni, tjID);
        end
    end

    [xi_est, yi_est, ss_var, x_predOrig, y_predOrig]=getSScores1(res,  x, t_x, y, t_y, nsvd, sig1, error, method, shrink, regular, rho_opt);


end

% function [y_predOrig,xi_est]=getOriCurves(y, t, mu, phi, lambda, sigma, sig1, noeig, error, method, shrink, out1, regular)
% 
% [muSub, phiSub] = convertMuPhi(t, out1, mu, phi, regular);
% ncohort = length(y);
% LAMBDA = lambda;
% 
% if error==1
% 
%     sigma1 = sig1;
% 
%     if regular == 2 && strcmp(method,'CE')
% 
%         yy= reshape(cell2mat(y), length(y{1}), ncohort)';
% 
%         error0 = sigma1*eye(length(t{1}));
%         A = LAMBDA*phiSub'*pinv(phiSub*LAMBDA*phiSub'+error0);
%         MU = repmat(muSub, ncohort,1);
%         B = yy-MU;
%         xi_est = (A*B')';
%         y_predOrig = MU+xi_est*phiSub';
%         y_predOrig = num2cell(y_predOrig,2);
% 
%     else
%         y_predOrig = cell(1,ncohort);
%         xi_est = zeros(ncohort, noeig);
%         zeta_est = xi_est;
%         phii= phiSub;
%         mu_i = muSub;
%         for i = 1:ncohort
% 
%             if regular ~= 2
%                 phii = phiSub{i};
%                 mu_i = muSub{i};
%             end
%             yi= y{i};
%             if strcmp(method,'CE')
%                 error0=sigma1*eye(length(yi));
%                 A = LAMBDA*phii'*pinv(phii*LAMBDA*phii'+error0);
%                 xi_est(i,:)=(A*(yi-mu_i)')';
% 
%             elseif strcmp(method,'IN')
%                 m=length(yi);
%                 for k=1:noeig
%                     prod=(yi-mu_i).*phii(:,k)';
%                     if shrink == 0
%                         %xi_est(i,k) = romb(t{i},prod);
%                         xi_est(i,k) = trapz(t{i},prod);
%                     else
%                         %zeta_est(i,k) = romb(t{i},prod);
%                         zeta_est(i,k) = trapz(t{i},prod);
%                         xi_est(i,k)=lambda(k)*zeta_est(i,k)/(lambda(k)+sigma/m);
%                     end
%                 end
%             end
%             y_predOrig{i} = mu_i+xi_est(i,:)*phii';
%         end
%     end
% elseif error==0
%     if regular == 2 && strcmp(method ,'CE')
%         yy= reshape(cell2mat(y), length(y{1}), ncohort)';
%         A = LAMBDA*phiSub'*pinv(phiSub*LAMBDA*phiSub');
%         MU = repmat(muSub, ncohort,1);
%         B = yy-MU;
%         xi_est = (A*B')';
%         y_predOrig = MU+xi_est*phiSub';
%         y_predOrig = num2cell(y_predOrig,2);
%     else
%         y_predOrig = cell(1,ncohort);
%         xi_est = zeros(ncohort, noeig);
%         phii= phiSub;
%         mu_i = muSub;
%         for i = 1:ncohort
%             if regular ~= 2
%                 phii = phiSub{i};
%                 mu_i = muSub{i};
%             end
%             yi= y{i};
%             if strcmp(method,'CE')
%                 A = LAMBDA*phii'*pinv(phii*LAMBDA*phii');
%                 xi_est(i,:)=(A*(yi-mu_i)')';
%             elseif strcmp(method,'IN')
%                 for k=1:noeig
%                     prod=(yi-mu_i).*phii(:,k)';
%                     %xi_est(i,k) = romb(t{i},prod);
%                     xi_est(i,k) = trapz(t{i},prod);
%                 end
%             end
%             y_predOrig{i} = mu_i+xi_est(i,:)*phii';
%         end
%     end
% end
% end

