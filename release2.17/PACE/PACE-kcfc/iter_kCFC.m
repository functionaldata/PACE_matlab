function [reCluster,newCluster,idselect,SSE,groupM] = iter_kCFC(y,t,M,initCluster,iter,op_cv,moveprop);
% iter_KCFC
%
% Last modified: 2015/3/8 by Pai-Ling Li
%
    n = length(y);
    out1 = unique(cell2mat(t));
    uniqclust = unique(initCluster);
    nclust = length(uniqclust);
            
    % -------------------------------------------------------------------------  
    % Main Steps : Iterative updating and reclassification
    % -------------------------------------------------------------------------      
    reCluster = zeros(n,iter);
    SSE = zeros(iter+1,nclust);
    groupM = SSE;
    idc = initCluster;
    
    for kk = 1:iter
        disp(['iteration :',num2str(kk)])
        idold = idc;
        %display(find(initCluster ~= idold))
      
        % ----------------------------------------------------------------------
        % Step 1: Estimate the model components and calculate the SSE based on
        % the clustering result of the (kk-1)th iteration.
        % ----------------------------------------------------------------------
        [groupM(kk,:),groupMufcn,groupEigfcn,groupEigval,groupevar,SSE(kk,:)] = allsse(y,t,nclust,idold,uniqclust,M);

        % ----------------------------------------------------------------------
        % Step 2 :  Reclassify the n curves
        % ----------------------------------------------------------------------
        diff = 10000000 * ones(n,nclust);
        if op_cv == 0
            % -------------------------------------------------------------------
            % Step 2-1 :  Reclassification without cross-validation (op_cv == 0) 
            % -------------------------------------------------------------------
            display(['Reclustering without CV...']);
            Ypred = cell(1,nclust);
            for k = 1:nclust
                if length(find(idold == uniqclust(k))) <= 1
                    continue; 
                end
                npcs = groupM(kk,k);
                mufcn = groupMufcn{k};
                tempEigfcn = groupEigfcn{k};
                lambda = groupEigval{k};
                evar = groupevar(k);
                [Fpcs, ~, Ypred{k}] = getScores(y,t,mufcn,tempEigfcn,lambda,evar,npcs, ...
                                        M.error,M.method,M.shrink,out1,M.regular,M.rho,M.verbose);
                for i = 1:n
                    %diff(i,k) = sum((y{i}-Ypred{k}{i}).^2);
                    diff(i,k) = trapz(t{i},(y{i}-Ypred{k}{i}).^2);
                end
            end;   %  k           
           
            for i = 1:n;
                tempdiff = diff(i,:);
                idtemp = find( tempdiff == min(tempdiff) ); 
                idc(i) = idtemp(1);                          
            end;
        elseif op_cv == 1; 
            % -------------------------------------------------------------------
            % Step 2-2 :  Reclassification with cross-validation (op_cv == 1)    
            % -------------------------------------------------------------------
            display(['Reclustering with CV...']);
            M_temp = M; % used for reclassification since CV need a fixed number of components.
            for i = 1:n
                Ypred = cell(1,nclust);
                for k = 1:nclust
                    display(['Reclustering with CV... subj ' num2str(i) ' cid ' num2str(idold(i)) ' on clust ' num2str(k)]);
                    if length(find(idold == uniqclust(k))) <= 1 
                        continue; 
                    end
                    npcs = groupM(kk,k);
                    M_temp.selection_k = npcs;

                    if idold(i) == k    % If the ith curve belongs to the kth cluster.
                        id = find(idold == uniqclust(k));
                        id(id == i) = [];
                        ncurve = length(id);
                        if ncurve <= 1 
                            continue; 
                        end
                        y_temp = y(id);
                        t_temp = t(id); 
                        M_temp.newdata = out1;
                        %yy_temp = FPCA(y_temp, t_temp, M_temp);
                        yy_temp = FPCAw(y_temp, t_temp, M_temp);
                        mufcn = getVal(yy_temp, 'mu');
                        %covfcn = getVal(yy_temp, 'xcov');
                        tempEigfcn = getVal(yy_temp, 'phi');
                        tempEigval = getVal(yy_temp, 'lambda');
                        tempFpcs = getVal(yy_temp, 'xi_est');
                        %varprop = getVal(yy_temp, 'FVE');
                        evar = getVal(yy_temp, 'sigma');
                    else % If the ith curve does not belong to the kth cluster.
                        id = find(idold == uniqclust(k));
                        mufcn = groupMufcn{k};
                        tempEigfcn = groupEigfcn{k};
                        tempEigval = groupEigval{k};
                        evar = groupevar(k);
                    end % idold(i) == k    

                    [fpcs_i, ~, Ypred(k)] = getScores(y(i),t(i),mufcn,tempEigfcn,tempEigval,evar,npcs, ...
                                                M.error,M.method,M.shrink,out1,M.regular,M.rho,M.verbose);
                    %diff(i,k) = sum((y{i}-Ypred{k}).^2);
                    diff(i,k) = trapz(t{i},(y{i}-Ypred{k}).^2);
                end   %  k
                tempdiff = diff(i,:);
                idtemp = find( tempdiff == min(tempdiff) ); 
                idc(i) = idtemp(1);
            end   %  i
        end   % op_cv
    
        % ----------------------------------------------------------------------
        % Step 3 : Check if the current clustering result is repeated in
        % any previous iteration.
        % ----------------------------------------------------------------------
        reCluster(:,kk) = idc;
        newCluster = idc;
    
        true = 0;
        if initCluster == newCluster
            if kk == 1;   
                true = 1;
            else;         
                true = 2;
            end;
            disp('Warning - KCFC : repeated clustering result = initial');
        else
            for j = 1:kk-2
                if reCluster(:,j) == newCluster
                    true = 2;
                    disp(['Warning - KCFC : repeated clustering result = iteration ',num2str(j)]);
                    break
                end
            end
        end
        if true == 1 
            idselect = 1; 
            break; 
        elseif  true == 2
            idnan = find(SSE == -99);
            if length(idnan) > 0
                tempSSE = SSE; 
                tempSSE(idnan) = 0;
                sumSSE = sum(tempSSE,2);
            else
                sumSSE = sum(SSE,2); 
            end   
            idselect = find(sumSSE(1:kk-1) == min(sumSSE(1:kk-1)));
            idselect = idselect(1);
            if idselect == 1
                newCluster = initCluster;
            else  
                newCluster = reCluster(:,idselect-1);
            end
            return;
        end
     
        % ----------------------------------------------------------------------
        % Step 4 :  Check if the stopping rule is satisfied.
        % ----------------------------------------------------------------------
        idmove = logical(idc == idold);
        nmove = length(find(idmove == 0));
        if nmove <= moveprop*n
            disp(['kCFC iteration converged after ', num2str(kk), ' reallocations.']);
            idselect = kk + 1;
            [groupM(kk+1,:),groupMufcn,groupEigfcn,groupEigval,groupevar,SSE(kk+1,:)] = allsse(y,t,nclust,idc,uniqclust,M);
            groupevar
            return;
        end
    
        if kk == iter
            idselect = kk + 1;
            [groupM(kk+1,:),groupMufcn,groupEigfcn,groupEigval,groupevar,SSE(kk+1,:)] = allsse(y,t,nclust,idc,uniqclust,M);
            disp(['Warning - KCFC : no covergent result after ', num2str(iter), ' reallocations.']);
        end
    end % end of kk
end
