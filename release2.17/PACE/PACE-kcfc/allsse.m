function [groupM,groupMufcn,groupEigfcn,groupEigval,groupevar,groupSSE] = allsse(y,t,nclust,idclust,uniqclust,M)
% ALLSSE
%
% Last modified: 2015/3/5 by Pai-Ling Li
%
    %n = size(y,2);
    groupM = zeros(1,nclust);
    groupMufcn = cell(1,nclust);
    groupEigfcn = cell(1,nclust);
    groupEigval = cell(1,nclust);
    groupevar = zeros(1,nclust);
    groupSSE = zeros(1,nclust);

    for k = 1:nclust;
        display(['Computing SSE on cluster ', num2str(k)]);
        id = find(idclust == uniqclust(k));   
        ncurve = length(id);
        if ncurve <= 1;
            groupSSE(k) = -99;
            gruopM(k) = -99;
        elseif ncurve > 1;
            y_temp = y(id);
            t_temp = t(id);
            %yy_temp = FPCA(y_temp,t_temp,M);
            yy_temp = FPCAw(y_temp,t_temp,M);
            groupM(k) = getVal(yy_temp,'no_opt');
            groupMufcn{k} = getVal(yy_temp,'mu');
            groupEigfcn{k} = getVal(yy_temp,'phi');
            groupEigval{k} = getVal(yy_temp,'lambda');
            if M.error == 0;
               groupevar(k) = 0;
            elseif M.error == 1;
               groupevar(k) = getVal(yy_temp,'sigma');
            end;
            %covfcn = getVal(yy_temp,'xcov');
            %PCs = getVal(yy_temp,'xi_est');
            ypred_temp = getVal(yy_temp,'y_predOrig');
            tempsse = 0;
            for i = 1:ncurve
                %tempsse = tempsse + sum((y_temp{i}-ypred_temp{i}).^2);
                tempsse = tempsse + trapz(t_temp{i},(y_temp{i}-ypred_temp{i}).^2);
            end
            groupSSE(k) = tempsse;
       end;  % ncurve > 1
    end;  %  k
end
